"""This script is to demonstrate Maximum Likelihood Estimation"""
# Created: Feb 2024, C.R. Zhao

# ---------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import numpy                as np
import matplotlib.pyplot    as plt
from scipy.optimize         import minimize
from sklearn.preprocessing  import PolynomialFeatures
from sklearn.linear_model   import LinearRegression, HuberRegressor
from sklearn.pipeline       import Pipeline
from sklearn.metrics        import mean_squared_error, mean_absolute_error


# ---------------------------------------------------------------------
#   Functions
# ---------------------------------------------------------------------
def GetPolynomialData():
    # generate outputs of a defined polynomial
    x       = np.linspace(-4.0, 4.0, num=100)
    a       = -1.0
    b       = 0
    c       = 9.0
    d       = 0.0
    e       = 0.0
    g       = 0.0
    y       = (a*(x**5)) + (b*(x**4)) + (c*(x**3)) + (d*(x**2)) + (e*x) + g

    # add noise to the data
    # np.random.normal(mean, standardDeviation, num)
    noise   = np.random.normal(0.0, 30.0, 100)
    y       = y + noise

    # save the data as a np array
    np.save('PolynomialData', [x, y])


def plot_PolynomialData():
    [x, y]  = np.load('PolynomialData.npy')
    fig1    = plt.figure()
    plt.plot(x, y, 'ok')
    plt.grid(True)
    plt.show()
    fig1.savefig('MLE_data.png', dpi=300, format='png')


# define a function to calculate the log likelihood
def calcLogLikelihood(guess, true, n):
    error   = true-guess
    sigma   = np.std(error)
    f       = ((1.0 / (2.0 * np.pi * sigma * sigma))**(n / 2)) * \
              np.exp(-1 * ((np.dot(error.T, error))/(2 * sigma * sigma)))
    return np.log(f)


# define my function which will return the objective function to be minimized
def OptFunction(coeffs, degree):
    [x, y]  = np.load('PolynomialData.npy')
    yGuess  = 0
    for i in range(degree):
        yGuess  += coeffs[i] * x**(i)
    f       = calcLogLikelihood(yGuess, y, float(len(yGuess)))
    return (-1 * f)


# calculate the results and plot MLE regression
def runOptimisationMLE():
    # Pick some random starting points for the optimisation
    degree      = 6
    coeffs      = np.zeros(degree)
    coeffs[0]   = -12.0
    coeffs[1]   = 18.0
    coeffs[2]   = -1.5
    coeffs[3]   = 1.0
    coeffs[4]   = 2.0
    coeffs[5]   = -1.0
    # methods: BFGS, SLSQP, CG, Powell...
    res         = minimize(OptFunction, coeffs, args=(degree), method='BFGS', options={'disp': True})

    # data sets
    [x, y]      = np.load('PolynomialData.npy')

    # other regression models
    # Linear Regressor
    model1 = Pipeline([('poly', PolynomialFeatures(degree=6)), ('linear', LinearRegression(fit_intercept=False))])
    model1 = model1.fit(x[:, np.newaxis], y)

    # calculate error metrics
    y_pred_MLE      = (res.x[5]*(x**5)) + (res.x[4]*(x**4)) + (res.x[3]*(x**3)) + (res.x[2]*(x**2)) + (res.x[1]*x) + res.x[0]
    y_pred_LR       = model1.predict(x[:, np.newaxis])

    # MSE
    err_MSE_MLE     = mean_squared_error(y, y_pred_MLE)
    err_MSE_LR      = mean_squared_error(y, y_pred_LR)
    # RMSE
    err_RMSE_MLE    = mean_squared_error(y, y_pred_MLE, squared=False)
    err_RMSE_LR     = mean_squared_error(y, y_pred_LR, squared=False)
    # MSE
    err_MAE_MLE     = mean_absolute_error(y, y_pred_MLE)
    err_MAE_LR      = mean_absolute_error(y, y_pred_LR)
    print('**********MLE method***************')
    print('err_MSE_MLE = {}'.format(err_MSE_MLE), '; err_RMSE_MLE = {} '.format(err_RMSE_MLE), '; err_MAE_MLE = {} '.format(err_MAE_MLE))
    print('**********LR method***************')
    print('err_MSE_LR = {}'.format(err_MSE_LR), '; err_RMSE_LR = {} '.format(err_RMSE_LR), '; err_MAE_LR = {} '.format(err_MAE_LR))

    # plot the data and model fits
    fig2    = plt.figure()
    plt.plot(x, y, 'ok')
    plt.plot(x, y_pred_MLE, '--b', label='MLE method')
    plt.plot(x, y_pred_LR, '.r', label='LR method')
    plt.grid(True)
    plt.legend(loc='best')
    plt.show()
    fig2.savefig('MLE_fitting.png', dpi=300, format='png')


# ---------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------
if __name__ == '__main__':
    # step 1: generate results
    GetPolynomialData()

    # step 2: plot data with noise
    plot_PolynomialData()

    # step 3: inversely solve the problem to determine the parameters using MLE method and plot it
    runOptimisationMLE()

