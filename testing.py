from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
import seaborn as sns
import pandas_datareader.data as web
import pandas as pd
from datetime import datetime
import scipy.stats as stats

c_funcs = CDLL("c_assets/cmake-build-debug/libMonteCarlo.dylib")


def c_geometricBrownianMT(n_steps, n_sims, max_t, W0, S0, mu,sigma,n_threads=1):
    _gmbMT = c_funcs._gmbThreadedInit
    _gmbMT.argtypes = [c_int, c_int, c_int, c_double, c_double, c_double, c_double, ndpointer(dtype=c_double, ndim=1, shape=(n_sims*n_steps))]
    _gmbMT.restype = None
    ret = np.empty(n_steps*n_sims, dtype=c_double)
    _gmbMT(n_sims, n_steps, n_threads, S0, W0, mu, sigma, ret)
    return ret


def extractParams(stock_name, start_date, end_date):
    df = web.DataReader(stock_name, data_source="yahoo", start=start_date, end=end_date)
    adj_close = df["Adj Close"]
    daily_returns = np.array([(y - x) / x for x, y in zip(adj_close, adj_close[1:])])
    n_days = len(df)
    mu = daily_returns.mean() * n_days
    sigma = daily_returns.std() * np.sqrt(n_days)
    return mu, sigma, n_days


def main():
    stock = 'GOOG'
    sns.set(style="darkgrid")
    n_sim = 1000
    mu, sigma, days = extractParams(stock, datetime(2019,1,31), datetime(2019,12,31))
    df = web.DataReader(stock, data_source="yahoo", start=datetime(2020,1,1), end=datetime.now())
    adj_close = df["Adj Close"]
    daily_returns = np.array([(y - x) / x for x, y in zip(adj_close, adj_close[1:])])
    n_steps = len(adj_close)
    S = c_geometricBrownianMT(n_steps, n_sim, 1, 0, adj_close[0], mu, sigma, n_threads=12)
    S = np.reshape(S, (-1, n_steps))
    fig, ax = plt.subplots(2,2,dpi=300)
    for i in range(n_sim):
        if i % 50 == 0:
            ax[0,0].plot(S[i, :], lw=.5)
    ax[0,0].plot(adj_close.values, lw=1.5, color="red", label='actual stock')
    ax[0,0].legend()
    ax[0,0].set_title('Simulated trajectories')
    ax[0,1].set_title('Stock price distribution for day {0}'.format(n_steps))
    ax[1,0].set_title('Trajectory average')
    sns.histplot(S[:, -1], kde=True, ax=ax[0,1])
    kde = stats.gaussian_kde(S[:, -1])
    price_range = np.linspace(S[:, -1].min(), S[:, -1].max(), n_sim)
    ax[0,1].axvline(price_range[np.argmax(kde(price_range))], c='green', ls = '-.', label='most likely price')
    ax[0,1].axvline(adj_close[-1], color='purple', ls='--', label = 'actual price')
    print("Prediction delta:", price_range[np.argmax(kde(price_range))]-adj_close[-1])
    print("relative prediction error:", 100*(price_range[np.argmax(kde(price_range))]-adj_close[-1])/adj_close[-1])
    ax[0,1].legend()
    traj_mean = S.mean(axis=0)
    traj_std = S.std(axis=0)
    ax[1,0].plot(traj_mean, c='k', ls='--', label='mean')
    ax[1,0].fill_between(range(n_steps), traj_mean+traj_std, traj_mean-traj_std, color='yellow', label='stand. deviation', alpha=0.4)
    ax[1, 0].plot(adj_close.values, lw=1, color="red")
    ax[1,0].legend()
    mlp = []
    for i in range(n_steps):
        kde = stats.gaussian_kde(S[:, i])
        price_range = np.linspace(S[:, i].min(), S[:, i].max(), n_sim)
        mlp.append(price_range[np.argmax(kde(price_range))])
    ax[1,1].plot(mlp, label="most likely trajectory")
    ax[1, 1].plot(adj_close.values, lw=1, color="red", label="actual stock")
    ax[1,1].set_title("Most likely trajectory")
    ax[1,1].legend()
    plt.show()


if __name__ == "__main__":
    main()