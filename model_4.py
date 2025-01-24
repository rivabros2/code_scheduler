import pandas as pd
import numpy as np
from sklearn import model_selection, ensemble, linear_model, neural_network, metrics, inspection

import matplotlib.pyplot as plt
from combat.pycombat import pycombat
from pickle import load, dump

# Para soyviz (contiene funciones creadas por Guido en su repositorio):
import sys
sys.path.append('/home/shared/epilab/BiomarkersML/soy_biomarkers_Guido/modelado/')

from soyviz import plot_los_tres, plot_estres

directorio_en_el_cluster = '/home/shared/epilab/BiomarkersML/'

venancio_genes = pd.read_parquet(directorio_en_el_cluster + 'SoyBean_MLresults/venancio_genes_bias_corrected.parquet')
print('venenacio genes:\n', venancio_genes)


annot = pd.read_csv(directorio_en_el_cluster + 'soy_biomarkers_Guido/anotacion/gemini_annot1.csv', index_col=0)
annot = annot[annot['tejido'].notna() & annot['estres'].notna()]
# annot = annot.loc[(annot['tejido'] == 'leaf') | (annot['tejido'] == 'root') | (annot['tejido'] == 'seed') | (annot['tejido'] == 'stem')]
print('annot:\n', annot)


venancio_genes, annot = venancio_genes.align(annot, join='inner', axis=0)
print('venancio:\n', venancio_genes)


consamples = annot['bioproject'].map(annot['bioproject'].value_counts() > 2)
convarianza = venancio_genes.var() > 0.01

filtrado = venancio_genes.loc[consamples, convarianza]

print('filtrado:\n', filtrado)


filtrado, annot = filtrado.align(annot, join='inner', axis=0)

combat = pycombat(filtrado.T, annot['bioproject']).T
print('combat:\n', combat)



# ---------- Separacion de los datos ----------

X, y = combat.align(annot, join='inner', axis=0)
print('X:\n', X)

plot_los_tres(X, y)
plt.suptitle('Muestras de PRJNA706999')
plt.savefig(directorio_en_el_cluster + 'SoyBean_MLresults/graficos/plot_estres_genes_04.png')
#plt.show()

y = y['estres']
print('conteo de valores de y:\n', y.value_counts())


X_train, X_test, y_train, y_test = model_selection.train_test_split(X, y)
print('X_train:\n', X_train)

# ---------- RANDOM FOREST (con Cross Validation) ----------

random_forest = ensemble.RandomForestClassifier()
params = {
    'n_estimators': np.arange(100, 1000, 100),
    'max_depth': np.arange(20, 50, 10),
    'min_samples_split': np.arange(2, 15, 4),
    'min_samples_leaf': np.arange(2, 15, 4),
    'max_features': ['sqrt', 'log2'],
    'bootstrap': [True]
}

from sklearn.model_selection import RandomizedSearchCV

try:
	with open('SoyBean_MLresults/forest_4.pkl', 'rb') as f:
		forest = load(f)
		print("si")
except FileNotFoundError:
	search = model_selection.RandomizedSearchCV(random_forest, param_distributions=params, n_jobs=-1, scoring='roc_auc_ovr', verbose=3)
	search.fit(X_train, y_train)
	
	forest = search.best_estimator_
	print("no")
	with open('SoyBean_MLresults/forest_4.pkl', 'wb') as f:
		dump(forest, f, protocol=5)


forest = search.best_estimator_
forest.get_params()

forest.score(X_test, y_test)

metrics.f1_score(y_test, forest.predict(X_test), average='macro')


importances = pd.DataFrame({'importance': np.flip(np.sort(forest.feature_importances_))}, index=X_train.columns[np.flip(forest.feature_importances_.argsort())])
plt.semilogx(importances['importance'].to_list(), '.')
plt.title('Importancias de features ordenadas (Random Forest)')
plt.savefig(directorio_en_el_cluster + 'SoyBean_MLresults/graficos/feature_importances_genes_04.png')
#plt.show()

importances.to_csv(directorio_en_el_cluster + 'SoyBean_MLresults/importancias_04_genes_rf.csv')


from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import seaborn as sns

y_pred = forest.predict(X_test)

cm = confusion_matrix(y_test, y_pred)

labels = list(np.unique(y_test))

# Configura la visualización de la matriz de confusión
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap="crest", vmin=0, linewidth=.5, xticklabels=labels, yticklabels=labels)
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion Matrix')
plt.savefig(directorio_en_el_cluster + 'SoyBean_MLresults/confusionmatrix_04_rf_genes.png')
##plt.show()




categories = y_test.unique()
label_mapping = {category: idx for idx, category in enumerate(categories)}
y_test_mapped = y_test.map(label_mapping)
y_pred_class = pd.Series(y_pred).map(label_mapping)



# Generar el reporte de clasificación
print('Reporte de la clasificación:\n', metrics.classification_report(y_test_mapped, y_pred_class, target_names=categories))



# ---------- GRADIENT BOOSTING (con Cross Validation) ----------


gbr = ensemble.GradientBoostingClassifier()
params = {
    "loss": ['log_loss'],
    "learning_rate": np.arange(0, 10, 1),
    "n_estimators": np.arange(1, 200, 1),
    "max_depth": np.arange(1, 50, 1),
    "max_features": ['sqrt', 'log2']
}


#gbr = ensemble.GradientBoostingClassifier()
#params_gbr = {
#    "loss": ['log_loss'],
#    "learning_rate": np.linspace(0.01, 0.2, 10),  # Rango más chico para tasas de aprendizaje
#    "n_estimators": np.arange(50, 300, 50),       
#    "max_depth": np.arange(3, 10, 1),             # Profundidad limitada no sobreajustar
#    "max_features": ['sqrt', 'log2']
#}

try:
	with open('SoyBean_MLresults/gradient_4.pkl', 'rb') as f:
		gradient = load(f)
		print("si")
except FileNotFoundError:
		search = model_selection.RandomizedSearchCV(gbr, param_distributions=params, n_jobs=-1, scoring='roc_auc_ovr', verbose=3)
		search.fit(X_train, y_train)
		gradient = search.best_estimator_
		print("no")
		with open('SoyBean_MLresults/gradient_4.pkl', 'wb') as f:
			dump(gradient, f, protocol=5)






gradient.score(X_test, y_test)
metrics.f1_score(y_test, gradient.predict(X_test), average='macro')



importances_gbr = pd.DataFrame({'importance': np.flip(np.sort(gradient.feature_importances_))}, index=X_train.columns[np.flip(gradient.feature_importances_.argsort())])
importances_gbr.to_csv(directorio_en_el_cluster + 'SoyBean_MLresults/importancias_04_gbr_genes.csv')

plt.semilogx(importances_gbr['importance'].to_list(), '.')
plt.title('Importancias de features ordenadas (Gradient boosted trees)')
plt.savefig(directorio_en_el_cluster + 'SoyBean_MLresults/feature_importances_04_gbr_genes.png')
#plt.show()



import seaborn as sns
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

y_pred = gradient.predict(X_test)

cm = confusion_matrix(y_test, y_pred)

labels = list(np.unique(y_test))

# Configura la visualización de la matriz de confusión
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap="crest", vmin=0, linewidth=.5, xticklabels=labels, yticklabels=labels)
plt.xlabel('Predicted')
plt.ylabel('True')
plt.xticks(fontsize=8, rotation=45, ha="right")
plt.yticks(fontsize=8)
plt.title('Confusion Matrix')
plt.savefig(directorio_en_el_cluster + 'SoyBean_MLresults/confusionmatrix_04_grb_genes.png')
#plt.show()

categories = y_test.unique()
label_mapping = {category: idx for idx, category in enumerate(categories)}
y_test_mapped = y_test.map(label_mapping)
y_pred_class = pd.Series(y_pred).map(label_mapping)

# Generar el reporte de clasificación
print('Reporte de la clasificación:\n', metrics.classification_report(y_test_mapped, y_pred_class, target_names=categories))


# (...)
