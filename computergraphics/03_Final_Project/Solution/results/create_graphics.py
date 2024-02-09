import json
import pandas as pd
import matplotlib.pyplot as plt

times = pd.read_json('data/render.times')
psnr = pd.read_json('data/render.psnr')

times.plot.line(grid=True, yticks=[0, 100, 200, 300, 400, 500], xlabel='samples', ylabel='time(s)', title='Time to render')
plt.savefig('graphs/times.png')
times.plot.line(subplots=True, grid=True, yticks=[0, 250, 500], ylim=(0, 600), xlabel='samples', ylabel='time(s)', title='Time to render')
plt.legend(loc='upper left')
plt.savefig('graphs/times_split.png')


psnr.plot.line(grid=True, ylim=(30, 55), xlabel='samples', ylabel='PSNR', title='Peak Signal to Noise ratio (higher=better)')
plt.legend(loc='lower right')
plt.savefig('graphs/psnr.png')
psnr.plot.line(subplots=True, grid=True, ylim=(30, 55), xlabel='samples', ylabel='PSNR', title='Peak Signal to Noise ratio (higher=better)')
plt.savefig('graphs/psnr_split.png')