# coding=utf-8

import sys
import numpy as np
import uncertainties as unc
from uncertainties.unumpy import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


gamma = 0.24


def loadData(p1, p2):
    BetaX = []
    with open(p1, 'r') as read_f:
        for row in read_f:
            temp = []
            row = row.strip('\n').split('\t')
            for e in row:
                temp.append(float(e))
            BetaX.append(temp)

    truncation = []
    with open(p2, 'r') as read_f:
        for row in read_f:
            row = row.strip('\n')
            truncation.append(int(row))
    return BetaX, truncation


def my_fun(x, alpha, omega):
    """BetaX equation with two parameters (alpha, omega) when fixed gamma"""
    y = alpha * x * (np.power((1.0 - gamma), np.power(x, omega)))
    return y


def main():
    p1 = './BetaX.csv'
    p2 = './truncation.csv'
    BetaX, truncation = loadData(p1, p2)

    for mid in range(1, 11):
        x = [float(i) for i in range(1, truncation[mid - 1] + 1)]
        x = np.array(x)
        y = BetaX[mid - 1][0:len(x)]
        y = np.array(y)

        p_est, err_est = curve_fit(my_fun, x, y)

        alpha_est, omega_est = unc.correlated_values(p_est, err_est)
        y_est = alpha_est * x * (np.power((1.0 - 0.24), np.power(x, omega_est)))
        nom = unc.unumpy.nominal_values(y_est)
        std = unc.unumpy.std_devs(y_est)
        lb, ub = nom - 2 * std, nom + 2 * std

        plt.plot(x, y, "bo", label='Empirical')
        plt.plot(x, nom, "r-", label='Fitting')
        plt.fill_between(x, lb, ub, color="k", alpha=0.05)

        plt.title(r'subfigure-'+str(mid), fontdict={'family': 'Times New Roman', 'size': 16})
        plt.xlabel(r'Number of Exposures, $x$', fontdict={'family': 'Times New Roman', 'size': 14})
        plt.ylabel(r'Retweeting Probability, $\beta(x)$', fontdict={'family': 'Times New Roman', 'size': 14})
        plt.legend(prop={'family': 'Times New Roman', 'size': 14})
        plt.text(2, 1.3*max(y), r'$\alpha=$'+str(round(alpha_est.nominal_value, 4)), fontsize=12, color='black')
        plt.text(2, 1.2 * max(y), r'$\omega=$'+str(round(omega_est.nominal_value, 2)), fontsize=12, color='black')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        plt.xlim(0, len(x) + 1)
        plt.ylim(0, 1.5 * np.max(y))
        plt.show()


if __name__ == '__main__':
    main()
