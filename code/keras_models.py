import pandas as pd
import numpy as np
import pickle
import tensorflow.keras as keras

from keras.callbacks import ReduceLROnPlateau, ModelCheckpoint
from keras.callbacks import TensorBoard
from keras.models import Model
from keras.layers import Input, Dense, BatchNormalization
from keras.callbacks import EarlyStopping
from model_routines import (f1, new_bac)

from sklearn.model_selection import train_test_split
from datetime import date


def encoding_mlp_generator(hidden_layers, desc_shape):
    bnorm2 = BatchNormalization(trainable=True)
    model_input = Input(shape=desc_shape)
    connect_to = model_input
    hidden_layer = Dense(hidden_layers,
                         activation="relu")(connect_to)
    connect_to = bnorm2(hidden_layer)
    output_t = Dense(1, activation="sigmoid",
                     name="output_t")(connect_to)

    model = Model(inputs=model_input,
                  outputs=output_t)
    model.compile(optimizer="adam",
                  loss="binary_crossentropy",
                  metrics=[f1, new_bac, tensorflow.keras.metrics.AUC(name="auc")],
                  run_eagerly=True)
    model.summary()

    return model


def new_generator(validation_set, batch_size, desc_dict):
    samples_per_epoch = np.array(validation_set).shape[0]
    number_of_batches = samples_per_epoch / batch_size
    counter = 0
    while True:
        x_data = validation_set[batch_size * counter:batch_size * (counter + 1)]
        x_batch = np.stack([desc_dict[x][0][0] for x in x_data])
        y_batch = np.stack([desc_dict[x][1] for x in x_data])
        counter += 1
        yield x_batch, y_batch
        if counter >= number_of_batches:
            counter = 0

def train_filter_keras_model(config) -> None:
    with open("1,986,447_hashed_fp_4096bitLength_2-4R_4nbp.pickle", "rb") as pkl:
        desc_dict = pickle.load(pkl)
        len_desc = len(desc_dict)
        print(len_desc)

    desc_shape = len(pd.DataFrame(desc_dict["1"][0]).iloc[0])
    indices = list(desc_dict.keys())

    callback = EarlyStopping(monitor="val_new_bac",
                             patience=3,
                             mode="max",
                             min_delta=0.0001)

    NAME = "Model_{}_descriptors_{}".format(len_desc, date.today())

    tensorboard_callback = TensorBoard(log_dir='logs/{}'.format(NAME), histogram_freq=1)
    sess = tensorflow.compat.v1.Session()
    # file_writer = tensorflow.summary.FileWriter('/my_log_dir_2/', sess.graph)

    mcp = ModelCheckpoint("model.h5",
                          save_best_only=True,
                          monitor="val_new_bac",
                          mode="max")

    rlr = ReduceLROnPlateau(monitor="val_new_bac",
                            factor=0.5,
                            patience=1,
                            min_lr=0.000001,
                            mode="max",
                            min_delta=0.0001,
                            verbose=1)

    pretrain, test = train_test_split(indices,
                                      test_size=0.1,
                                      random_state=42,
                                      shuffle=True)
    train, validation_set = train_test_split(pretrain,
                                             test_size=0.2,
                                             random_state=41)

    model = encoding_mlp_generator(1000, desc_shape)

    history = model.fit(new_generator(train, 200, desc_dict),
                        steps_per_epoch=len(train) // 200,
                        validation_data=new_generator(validation_set, 200, desc_dict),
                        validation_steps=len(validation_set) // 200,
                        callbacks=[callback, rlr, mcp, tensorboard_callback],
                        verbose=1,
                        epochs=100,
                        use_multiprocessing=False,
                        workers=1)

    prediction = model.evaluate(new_generator(test, 200, desc_dict),
                                batch_size=200,
                                steps=len(test) // 200,
                                verbose=1,
                                use_multiprocessing=False,
                                workers=1)