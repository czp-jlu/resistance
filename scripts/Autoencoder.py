from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, BatchNormalization
from tensorflow.keras.optimizers import Adam
import os
import pandas as pd
import numpy as np

# working directory
os.chdir("C:/Users/sf8642/Desktop/meine_paper/Paper_2_Carola/revision/upload_data")

# working directory
os.chdir("C:/Users/$_here_goes_your_directory_$")

marker = pd.read_csv("data/Marker_matrix.csv", sep=";")
marker = marker/2
marker = marker.astype(np.uint32)

# splts  = pd.read_csv("data/assignment_sets.txt", sep=" ")
# idx = splts.iloc[:,0] == "T"
# train = marker.iloc[idx.array,:]
# val   = marker.iloc[~idx.array,:]


### Using 250 units for the intermediate layer
encoder = Sequential()
encoder.add(Dense(units=4000, activation='relu'))
encoder.add(BatchNormalization())
encoder.add(Dense(units=1000, activation='relu'))
encoder.add(BatchNormalization())
encoder.add(Dense(units=250, activation='relu'))

decoder = Sequential()
decoder.add(Dense(units=1000, activation='relu'))
decoder.add(BatchNormalization())
decoder.add(Dense(units=4000, activation='relu'))
decoder.add(BatchNormalization())
decoder.add(Dense(units=16667, activation='sigmoid'))

autoencoder = Sequential([encoder, decoder])
autoencoder.compile(loss="mse", optimizer=Adam(learning_rate=0.0001), metrics="binary_crossentropy")

hist = autoencoder.fit(marker, marker, epochs=100)

encoded_data = encoder.predict(marker)

res = pd.DataFrame(encoded_data)
res.to_csv("AE_result.csv")
autoencoder.save("models/AE_model")
encoder.save("models/ENC_model")
pd.DataFrame(hist.history).to_csv("fit_history_AE.csv")
