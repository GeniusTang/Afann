import numpy as np
from sklearn.neural_network import MLPRegressor
from sklearn.base import BaseEstimator, TransformerMixin

class padding_MLPR(BaseEstimator, TransformerMixin):
    def __init__(self, method, padding_ratio = 2, hidden_layer_sizes = 2000, seed = 42):
        self.padding_ratio = padding_ratio
        self.hidden_layer_sizes = hidden_layer_sizes
        self.rng = np.random.RandomState(seed)
        self.seed = seed
        self.model = MLPRegressor(hidden_layer_sizes=self.hidden_layer_sizes, random_state=self.seed)
        self.model.n_layers_ = 3
        self.model.n_outputs_ = 1
        self.model.out_activation_ = 'identity'
        self.model.coefs_ = np.load('model/{}_coefs.npy'.format(method))
        self.model.intercepts_ = np.load('model/{}_intercepts.npy'.format(method))
 
    def fit(self, X, y=None):
        
        padding = int(X.shape[0] * self.padding_ratio / 2)
        
        X_A1 = np.c_[np.linspace(0, 1, padding), np.ones(padding), np.ones(padding)]
        y_A1 = X_fake2[:,0]

        X_A2 = np.c_[np.zeros(padding), self.rng.rand(padding), self.rng.rand(padding)]
        y_A2 = np.zeros(padding)

        final_X = np.vstack([X, X_A1, X_A2])
        final_y = np.hstack([y, y_A1, y_A2])
        
        self.model.fit(final_X, final_y)
        return self

    def predict(self, X):
        X = np.array(X)
        return (self.model.predict(X) + self.model.predict(X[:,[0,2,1]])) / 2

    def score(self, X, y=None):
        X = np.array(X)
        return spearman_r(self.model.predict(X), y)

if __name__ == '__main__':
    d2shepp_model = padding_MLPR(method='d2shepp')
    print(d2shepp_model.predict([[0.4,1,1]]))
