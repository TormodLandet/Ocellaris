import numpy
from matplotlib import pyplot
from matplotlib.widgets import Slider


def read_line_probe_file(file_name):
    with open(file_name, 'rt') as f:
        # Read header
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()
        line4 = f.readline()
        f.readline()  # line 5 is not interesting

        # Parse header
        assert line1.startswith('# Ocellaris line probe of the')
        field_name = line1[30:-7].strip()
        xpos = numpy.array([float(w) for w in line2.split()[3:]], dtype=float)
        ypos = numpy.array([float(w) for w in line3.split()[3:]], dtype=float)
        zpos = numpy.array([float(w) for w in line4.split()[3:]], dtype=float)
        N = len(xpos)
        assert len(ypos) == len(zpos) == N

        # Read time steps
        time = []
        values = []

        for line in f:
            wds = line.split()
            if len(wds) != N + 1:
                break
            time.append(float(wds[0]))
            values.append(numpy.array([float(w) for w in wds[1:]], dtype=float))
        time = numpy.array(time, dtype=float)

        return field_name, xpos, ypos, zpos, time, values


def plot_line_probe_file(file_name):
    field_name, xpos, ypos, zpos, time, values = read_line_probe_file(file_name)

    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    pyplot.subplots_adjust(left=0.15, bottom=0.25)

    dx = xpos[-1] - xpos[0]
    dy = ypos[-1] - ypos[0]
    dz = zpos[-1] - zpos[0]
    d, pos, xlabel = dx, xpos, 'x'
    if dy > d:
        d, pos, xlabel = dy, ypos, 'y'
    if dz > d:
        pos, xlabel = zpos, 'z'

    line, = pyplot.plot(pos, values[-1])

    ax_slider = fig.add_axes([0.15, 0.1, 0.65, 0.03])
    slider = Slider(ax_slider, 'time', time[0], time[-1], valinit=time[-1])

    def update(val):
        i = numpy.argmin(abs(time - val))
        y = values[i]
        line.set_ydata(y)
        # ax.set_ylim(y.min(), y.max())
        fig.canvas.draw_idle()

    slider.on_changed(update)
    update(slider.val)

    ax.set_title(
        'Line probe of %s, avg(pos)=[%.2f, %.2f, %.2f]'
        % (field_name, xpos.mean(), ypos.mean(), zpos.mean())
    )
    ax.set_xlabel(xlabel)

    return slider  # avoid garbage collection of the slider making it non responsive


if __name__ == '__main__':
    import sys

    line_probe_file_name = sys.argv[1]
    slider = plot_line_probe_file(line_probe_file_name)
    pyplot.show()
