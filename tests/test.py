# -*- coding: utf-8 -*-
import cpp_spline_interp
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

test_results = pd.DataFrame(
	columns=['interp_method', 'extrap_method', 'x', 'y']
)

X = [1, 2, 3, 4, 5]
Y = [0, 1, 5, 9, 10]
for interp_method in ["linear", "cubic"]:
	for extrap_method in ["linear", "flat"]:

		spline_interp = cpp_spline_interp.SplineInterpolate(
			interp_method,
			extrap_method,
			X,
			Y,
		)
		x_vec = np.linspace(0, 6, 100)

		for idx in range(len(x_vec)):
			y = spline_interp.F(x_vec[idx])
			test_results = test_results.append(
				{
					'interp_method': interp_method,
					'extrap_method': extrap_method,
					'x': x_vec[idx],
					'y': y,
				},
				ignore_index=True,
			)


g = sns.FacetGrid(test_results, row="extrap_method", col="interp_method", margin_titles=True)

g.map(sns.lineplot, "x", "y", color=".3",)

g.axes[0, 0].plot(X, Y, 'ro')
g.axes[0, 1].plot(X, Y, 'ro')
g.axes[1, 0].plot(X, Y, 'ro')
g.axes[1, 1].plot(X, Y, 'ro')
plt.xticks([0, 1, 2, 3, 4, 5, 6])
plt.show()

