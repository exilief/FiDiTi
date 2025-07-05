import numpy as np
import matplotlib.pyplot as plot

outdir = "out"

default_width = plot.gcf().get_figwidth()

def plot_frame1D(i):
    file = outdir + "/fdtd1d." + str(i) + ".txt"
    data = np.fromfile(file, sep='\n')  # loadtxt
    xs = np.arange(0, len(data))

    #plot.gcf().set_figwidth(default_width / 1.5)
    #plot.gca().set_aspect(1)

    plot.title(r"$E_z$ (N/C)")
    plot.xlabel(r"$x / \Delta x$")
    #plot.ylabel(r"$E_z$ (N/C)")
    #plot.ylim(-data.max()*0.1, data.max()*1.1)
    plot.ylim(-0.1, 1.1)
    plot.tight_layout()
    plot.plot(xs, data) #(xs, data, 'm.')
    #plot.show()
    plot.savefig(outdir + "/plot/fdtd1d." + str(i) + ".png")

def plot_frame2D(i):
    file = outdir + "/fdtd2d." + str(i) + ".txt"
    data = np.loadtxt(file)
    m,n = len(data[0]), len(data)
    x,y = np.arange(0, m), np.arange(0, n)

    plot.title(r"$E_z$ (N/C)")
    plot.xlabel(r"$x / \Delta x$");  plot.ylabel(r"$y / \Delta y$")
    #plot.contourf(data)  # norm: linear, symlog, log
    if True:  # Fixed color range
        plot.contourf(np.round(data,5), 20, norm="symlog", vmin=-0.1, vmax=1) # (x, y, data, 16, norm="linear")
    else:
        plot.contourf(np.round(data,5), 20, norm="symlog")
    plot.colorbar()
    plot.savefig(outdir + "/plot/2d/fdtd." + str(i) + ".png")

def plot_frame2D_line(i):
    file = outdir + "/fdtd2d." + str(i) + ".txt"
    data = np.loadtxt(file)
    m,n = len(data[0]), len(data)
    x,y = np.arange(0, m), np.arange(0, n)

    plot.title('Ez (N/C)')
    plot.ylim(-0.1, 1.01)
    plot.plot(x, data[int(n/2)])
    plot.savefig(outdir + "/plot/2d/fdtd." + str(i) + ".png")

for i in range(57):
    plot_frame2D(i)
    plot.clf()


# Using np.abs(data) makes the TF/SF wave look much cleaner (no low-amplitude noise)
# (There are small oscillations around zero behind the wave -> Jumps between 2 colors)
# Alternative: np.round(data, 5)
