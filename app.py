# import pandas as pd
# df_pos = pd.read_csv('pos FVs.csv')

# df_neg = pd.read_csv('neg FVs.csv')
# df_neg = df_neg.sample(n=393, random_state=42)

# # assuming your two dataframes are called "df1" and "df2"
# df = pd.concat([df_pos, df_neg], ignore_index=True)

# import lightgbm as lgb
# from sklearn.metrics import accuracy_score, matthews_corrcoef, roc_curve, auc, confusion_matrix
# import joblib
# import numpy as np
# from sklearn.model_selection import train_test_split
# import pandas as pd
# from sklearn.preprocessing import StandardScaler
# from sklearn.preprocessing import MinMaxScaler

# # Assuming df is your DataFrame containing the data

# X = df.drop(['class'], axis=1).values
# y = df['class'].values

# std_scale = StandardScaler().fit(X)
# X = std_scale.transform(X)  # Only transform for consistency with test data
# X = np.nan_to_num(X.astype('float32'))

# # Split the dataset into training and testing sets
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

# # Define the LightGBM model architecture
# model = lgb.LGBMClassifier(random_state=42)

# # Train the model
# model.fit(X_train, y_train)

# # Save the trained model
# joblib.dump(model, 'lgbm_model.pkl')

# # Save the StandardScaler
# joblib.dump(std_scale, 'standard_scaler.pkl')

# # Predict on testing data
# y_pred = model.predict(X_test)

# # Calculate evaluation metrics
# accuracy = accuracy_score(y_test, y_pred)
# tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
# specificity = tn / (tn + fp)
# sensitivity = tp / (tp + fn)

# lgbm_fpr, lgbm_tpr, thresholds = roc_curve(y_test, y_pred)
# lgbm_roc_auc = auc(lgbm_fpr, lgbm_tpr)
# mcc = matthews_corrcoef(y_test, y_pred)

# # Print evaluation metrics
# print("Accuracy:", accuracy)
# print("MCC:", mcc)
# print("Specificity:", specificity)
# print("Sensitivity:", sensitivity)
# print("AUC-ROC:", lgbm_roc_auc)



# # Split sequences into individual entries
# entries = sequences.split('\n')

import re

entries = '''>sp|P10310|TERL_BPT3
MSTQSNRNALVVAQLKGDFVAFLFVLWKALNLPVPTKCQIDMAKVLANGDNKKFILQAFRGIGKSFITCAFVVWTLWRDPQLKILIVSASKERADLNSIFIKNIIDLLPFLDELKPSPGQRDSVISFDVGPAKPDHSPSVKSVGITGQLTGSRADIIIADDVEIPSNSATQGAREKLWTLVQEFRALLKPLPTSRVIYLGTPQTEMTLYKELEDNRGYTTIIWPALYPRSREEDLYYGERLAPMLREEFNDGFEMLQGQPTDPVRFDMEDLRERELEYGKAGFTLQFMLNPNLSDAEKYPLRLRDAIVCGLDFEKAPMHYQWLPNRQNRNEELPNVGLKGDDIHSYHSCSQNTGQYQQRILVIDPSGRGKDETGYAVLFTLNGYIYLMEAGGFPDGYSDKTLESLAKKANEWKVQTVVFESNFGDGMFGKVFSPVLLKHHAAALEEIRARGMKELRICDTLEPVLSTHRLVIRDEVIREDYQTARDADGKHDVRYSLFYQLTRMAREKGAVAHDDRLDAFRLGVEFLRSTMELDAVKVEAEVLEAFLEEHMEHPIHSAGHVVTAMVDGMELYWEDDDVNGDRFINW'''



# genename = entries.strip().split('|')[1].split('.')[0]


genename = entries.split('\n')[1]

print(genename)
# # Initialize a list to store the protein sequences
# protein_sequences = []

# # Iterate over each entry
# for entry in entries:
#     # Use regular expression to extract the protein sequence
#     matches = re.findall(r'[A-Z]+', entry)
#     protein_sequence = ''.join(matches[1:])  # Exclude the identifier
#     protein_sequences.append(protein_sequence)


# print(protein_sequences)