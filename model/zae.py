from keras.layers import Input, Dense, Dropout, Activation, BatchNormalization, Lambda
from keras.models import Model
from keras.regularizers import l1_l2
import keras.optimizers as opt
from keras.callbacks import ReduceLROnPlateau, EarlyStopping
from keras.objectives import mean_squared_error
import tensorflow as tf

from model import layers as custom_layers
from model import loss as custom_loss

MeanAct = lambda x: tf.clip_by_value(tf.math.exp(x), 1e-5, 1e6)
DispAct = lambda x: tf.clip_by_value(tf.nn.softplus(x), 1e-4, 1e4)

class ZINBAutoencoder():
    def __init__(self,
                 input_size,
                 output_size,
                 hidden_size=(64, 32, 64),
                 hidden_dropout=0.,
                 batchnorm=True,
                 activation='relu',
                 init='glorot_uniform',
                 debug=False):
      
        self.input_size = input_size
        self.output_size = output_size
        self.hidden_size = hidden_size
        self.hidden_dropout = hidden_dropout
        self.batchnorm = batchnorm
        self.activation = activation
        self.init = init
        self.debug = debug

        # potential args
        self.l2_coef = 0.
        self.l1_coef = 0.
        self.l2_enc_coef = 0.
        self.l1_enc_coef = 0.
        self.ridge = 0.
        self.input_dropout = 0.
        self.file_path = None

        # init
        self.extra_models = {}
        self.loss = None
        self.model = None
        self.encoder = None
        self.decoder = None
        self.input_layer = None
        self.sf_layer = None
        
        if isinstance(self.hidden_dropout, list):
            assert len(self.hidden_dropout) == len(self.hidden_size)
        else:
            self.hidden_dropout = [self.hidden_dropout]*len(self.hidden_size)

    def build(self):
        self.input_layer = Input(shape=(self.input_size,), name='count')
        self.sf_layer = Input(shape=(1,), name='size_factors')

        self.encoder = Dense(self.hidden_size[0], 
                             activation=None, 
                             kernel_initializer=self.init,
                             kernel_regularizer=l1_l2(self.l1_coef, self.l2_coef),
                             name="encoder")(self.input_layer)
        self.encoder = BatchNormalization(center=True, scale=False)(self.encoder)
        self.encoder = Activation(self.activation, name='act_encoder')(self.encoder)

        self.bottleneck = Dense(self.hidden_size[1], 
                             activation=None, 
                             kernel_initializer=self.init,
                             kernel_regularizer=l1_l2(self.l1_enc_coef, self.l2_enc_coef),
                             name="bottleneck")(self.encoder)
        self.bottleneck = BatchNormalization(center=True, scale=False)(self.bottleneck)
        self.bottleneck = Activation(self.activation, name='act_bottleneck')(self.bottleneck)

        self.decoder = Dense(self.hidden_size[2], 
                             activation=None, 
                             kernel_initializer=self.init,
                             kernel_regularizer=l1_l2(self.l1_coef, self.l2_coef),
                             name="decoder")(self.bottleneck)
        self.decoder = BatchNormalization(center=True, scale=False)(self.decoder)
        self.decoder = Activation(self.activation, name='act_decoder')(self.decoder)

        self.decoder_output = self.decoder
        self.build_output()

    def build_output(self):
        pi = Dense(self.output_size, activation='sigmoid', kernel_initializer=self.init,
                       kernel_regularizer=l1_l2(self.l1_coef, self.l2_coef),
                       name='pi')(self.decoder_output)

        disp = Dense(self.output_size, activation=DispAct,
                           kernel_initializer=self.init,
                           kernel_regularizer=l1_l2(self.l1_coef, self.l2_coef),
                           name='dispersion')(self.decoder_output)

        mean = Dense(self.output_size, activation=MeanAct, kernel_initializer=self.init,
                       kernel_regularizer=l1_l2(self.l1_coef, self.l2_coef),
                       name='mean')(self.decoder_output)
        output = custom_layers.ColwiseMultLayer([mean, self.sf_layer])
        output = custom_layers.SliceLayer(0, name='slice')([output, disp, pi])

        zinb = custom_loss.ZINB(pi, theta=disp, ridge_lambda=self.ridge, debug=self.debug)
        self.loss = zinb.loss
        self.extra_models['pi'] = Model(inputs=self.input_layer, outputs=pi)
        self.extra_models['dispersion'] = Model(inputs=self.input_layer, outputs=disp)
        self.extra_models['mean_norm'] = Model(inputs=self.input_layer, outputs=mean)
        self.extra_models['decoded'] = Model(inputs=self.input_layer, outputs=self.decoder_output)

        self.model = Model(inputs=[self.input_layer, self.sf_layer], outputs=output)

    def predict(self, adata, mode='denoise', return_info=False, copy=False, colnames=None):

        adata = adata.copy() if copy else adata

        if return_info:
            adata.obsm['X_dca_dispersion'] = self.extra_models['dispersion'].predict(adata.X)
            adata.obsm['X_dca_dropout']    = self.extra_models['pi'].predict(adata.X)

        print('dca: Calculating reconstructions...')

        adata.X = self.model.predict({'count': adata.X,
                                      'size_factors': adata.obs.size_factors})

        adata.uns['dca_loss'] = self.model.test_on_batch({'count': adata.X,
                                                          'size_factors': adata.obs.size_factors},
                                                          adata.raw.X)
        
        return adata if copy else None
                    
    def save(self):
        if self.file_path:
            os.makedirs(self.file_path, exist_ok=True)
            with open(os.path.join(self.file_path, 'model.pickle'), 'wb') as f:
                pickle.dump(self, f)

    def train(self, adata, optimizer='rmsprop', learning_rate=None,
          epochs=300, reduce_lr=10, use_raw_as_output=True,
          early_stop=15, batch_size=32, clip_grad=5., 
          validation_split=0.1, verbose=True, threads=None,
          **kwds):
        optimizer = opt.RMSprop(clipvalue=clip_grad)

        self.model.compile(loss=self.loss, optimizer=optimizer)

        # Callbacks
        callbacks = []

        if reduce_lr:
            lr_cb = ReduceLROnPlateau(monitor='val_loss', patience=reduce_lr, verbose=verbose)
            callbacks.append(lr_cb)
        if early_stop:
            es_cb = EarlyStopping(monitor='val_loss', patience=early_stop, verbose=verbose)
            callbacks.append(es_cb)

        if verbose: model.summary()

        inputs = {'count': adata.X, 'size_factors': adata.obs.size_factors}

        output = adata.raw.X if use_raw_as_output else adata.X

        loss = self.model.fit(inputs, output,
                        epochs=epochs,
                        batch_size=batch_size,
                        shuffle=True,
                        callbacks=callbacks,
                        validation_split=validation_split,
                        verbose=verbose,
                        **kwds)

        return loss