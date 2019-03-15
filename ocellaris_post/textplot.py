# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0
"""
Command line text plotting
==========================

Command line plotting tool inspired by hipsterplot, but with
significantly simplified code due to using numpy for histogram
construction.

The plot() function produces text plots like seen below from a list
of x and y coordinates. Scatter plots will also work to some degree,
as long as they are rather sparse::

    0.888                                        :.::              :
    0.584                      ::.:             .    .            :
    0.382                    ..    :           ..    ..          .
    0.179   ..::.:.         ..      :         .        .         .
  -0.0237 ::.     .:       ..        .        .        .        .
   -0.226           .:    :           .      .          .      .
   -0.429             :.::            ..    ..          .      .
   -0.632                              ..  :             ..   :
   -0.936                                ::               .::.
                          x axis range from 0 to 20

This is the result of the following code::

    x = numpy.linspace(0, 20, 100)
    y = numpy.sin(x) * numpy.tanh(x/10)
    plot(x, y, figsize=(70, 10))

The total length in characters including the y-axis ticks is here 70
and the total number of lines is 10 (including the line containing
the x-axis extents).
"""
import numpy


SYMBOLS = ' .:!||||++++'
INFINITY = '#'


def count_2_char(v):
    return SYMBOLS[v] if v < len(SYMBOLS) else INFINITY


def translate_line(vals):
    "Translate list of histogram counts to string representing density"
    return ''.join(count_2_char(v) for v in vals)


def plot(*args, figsize=(80, 15), xmin=None, xmax=None, ymin=None, ymax=None):
    """
    Plot a line or scatter plot

    Usage::

        plot(y)
        plot(x, y)

    You can specify width and height as size in characters to
    control the figure size. The with and height includes X and Y
    axis ticks, but no title or axis lables are shown
    """
    width, height = figsize
    assert width > 15
    assert height > 2

    # Unpack arguments
    if len(args) == 1:
        y = args[0]
        x = numpy.arange(len(y))
    elif len(args) == 2:
        x = args[0]
        y = args[1]
    else:
        # Fixme: we could support multiple lines with different ANSI colours
        raise NotImplementedError('You must give one or two positional arguments')

    # Filter out unwanted points (axes limits)
    if xmin is None:
        xmin = numpy.min(x)
    if xmax is None:
        xmax = numpy.max(x)
    if ymin is None:
        ymin = numpy.min(y)
    if ymax is None:
        ymax = numpy.max(y)
    x, y = filter_data_points(x, y, xmin, xmax, ymin, ymax)

    # Create 2D histogram of point density in each char "pixel"
    assert len(x) == len(y)
    W = width - 12
    H = height - 1
    hist, xedges, yedges = numpy.histogram2d(x, y, bins=(W, H))
    hist = numpy.asarray(hist, dtype=int)

    # Print 2D histogram representation of the line/scatter plot
    for i in list(range(H))[::-1]:
        if i == 0:
            ylabel = yedges[0]
        elif i == H - 1:
            ylabel = yedges[-1]
        else:
            ylabel = (yedges[i] + yedges[i + 1]) / 2
        print('% 9.3g' % ylabel, translate_line(hist[:, i]))
    xticks = 'x axis range from %g to %g (%d values)' % (xedges[0], xedges[-1], len(x))
    print(' ' * 9, xticks.center(W))


def filter_data_points(x, y, xmin, xmax, ymin, ymax):
    """
    Remove points outside the axes limits
    """
    xnew, ynew = [], []
    for xi, yi in zip(x, y):
        if xi is None or yi is None:
            continue
        elif xmin <= xi <= xmax and ymin <= yi <= ymax:
            xnew.append(xi)
            ynew.append(yi)
    return xnew, ynew
