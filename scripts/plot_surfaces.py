from __future__ import division, print_function
import sys
import numpy
from matplotlib import pyplot
from matplotlib.widgets import Slider


def read_iso_surfaces(filename):
    timesteps = []
    data = []

    with open(filename, 'rt') as f:
        description = f.readline()[1:].strip()
        value = float(f.readline().split()[-1])
        dim = int(f.readline().split()[-1])

        line = f.readline()
        while line:
            wds = line.split()
            try:
                time = float(wds[1])
                nsurf = int(wds[3])
            except Exception:
                break

            if nsurf == 0:
                timesteps.append(time)
                data.append([])
                line = f.readline()
                continue

            datalines = [f.readline() for _ in range(nsurf * 3)]
            if not datalines[-1]:
                break
            timesteps.append(time)
            data.append([])
            for i in range(nsurf):
                xvals = [float(v) for v in datalines[i * 3 + 0].split()]
                yvals = [float(v) for v in datalines[i * 3 + 1].split()]
                zvals = [float(v) for v in datalines[i * 3 + 2].split()]
                data[-1].append((xvals, yvals, zvals))

            line = f.readline()
    return description, value, dim, timesteps, data


def plotit(ax, contours, label):
    colour = get_colour(label)
    xvals = []
    yvals = []
    for contour in contours:
        if xvals:
            xvals.append(None)
            yvals.append(None)
        xvals.extend(contour[0])
        yvals.extend(contour[1])
    ax.plot(xvals, yvals, c=colour, label=label)


COLOURS = {}


def get_colour(label):
    available = 'bgrcmy'
    if label not in COLOURS:
        COLOURS[label] = available[len(COLOURS)]
    return COLOURS[label]


def main(filenames, equal_axes):
    all_data = []
    tmin = tmax = 0
    xmin, ymin, xmax, ymax = 1e100, 1e100, -1e100, -1e100
    for i, filename in enumerate(filenames):
        name = str(i + 1)
        if ':' in filename:
            filename, name = filename.split(':')
        print('Reading %s from file %s' % (name, filename))
        _description, _value, _dim, timesteps, data = read_iso_surfaces(filename)
        all_data.append((name, numpy.array(timesteps), data))

        tmax = max(tmax, timesteps[-1])
        for contours in data:
            for contour in contours:
                xmin = min(xmin, numpy.min(contour[0]))
                ymin = min(ymin, numpy.min(contour[1]))
                xmax = max(xmax, numpy.max(contour[0]))
                ymax = max(ymax, numpy.max(contour[1]))

    fig, ax = pyplot.subplots()
    pyplot.subplots_adjust(bottom=0.25)
    axcolor = '#a1b8dd'
    slider_ax = pyplot.axes([0.1, 0.1, 0.8, 0.03], facecolor=axcolor)
    slider = Slider(slider_ax, 'Time', tmin, tmax, valinit=tmin)

    xdiff = xmax - xmin
    ydiff = ymax - ymin
    xmin, xmax = xmin - 0.1 * xdiff, xmax + 0.1 * xdiff
    ymin, ymax = ymin - 0.1 * ydiff, ymax + 0.1 * ydiff
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    def update(val):
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.clear()

        t = slider.val
        for name, timesteps, data in all_data:
            i = numpy.argmin(abs(timesteps - t))
            dt = timesteps[1] - timesteps[0]
            if abs(timesteps[i] - t) > 1.5 * dt:
                continue
            plotit(ax, data[i], name)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        if equal_axes:
            ax.set_aspect('equal')

        if len(filenames) > 1:
            ax.legend(loc='lower right')
        fig.canvas.draw_idle()

    slider.on_changed(update)
    slider.set_val(tmin)
    pyplot.show()


if __name__ == '__main__':
    filenames = sys.argv[1:]

    equal_axes = True
    if '--nonequal' in filenames:
        equal_axes = False
        filenames.remove('--nonequal')

    main(filenames, equal_axes)
