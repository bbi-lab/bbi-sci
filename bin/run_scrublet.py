#!/usr/bin/env python
import matplotlib
matplotlib.use('pdf')
import scrublet as scr
import scipy.io
import numpy
import numpy.ma
from PIL import Image, ImageDraw, ImageFont
import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Check sample sheet')

    parser.add_argument('--key', required=True, help='Sample key (for naming).')
    parser.add_argument('--mat', required=True, help='Input matrix.')
    parser.add_argument('--skip', action="store_true", default=False, help='Should scrublet be skipped?')
    args = parser.parse_args()

    if not args.skip:
        counts_matrix = scipy.io.mmread(args.mat).T.tocsc()
        scrub = scr.Scrublet(counts_matrix)
        try:
            doublet_scores, predicted_doublets = scrub.scrub_doublets()
            scrub.plot_histogram()[0].savefig(args.key + "_scrublet_hist.png")
            all_scores = numpy.vstack((doublet_scores, predicted_doublets))
            all_scores = numpy.transpose(all_scores)
            numpy.savetxt(args.key + "_scrublet_out.csv", all_scores, delimiter=",")
        except (ZeroDivisionError, ValueError):
            temp = numpy.array(["NA"] * numpy.size(counts_matrix, 0))
            all_scores = numpy.vstack((temp, temp))
            all_scores = numpy.transpose(all_scores)
            filename = args.key + "_scrublet_hist.png"
            image = Image.new(mode = "RGB", size = (800,300), color = "white")
            draw = ImageDraw.Draw(image)
            draw.text((120,140), "Scrublet failed. This is generally because there aren't enough cells with sufficient reads.", fill = "black")
            image.save(filename)
            numpy.savetxt(args.key + "_scrublet_out.csv", all_scores, fmt="%s", delimiter=",")
        except (AttributeError):
            predicted_doublets = scrub.call_doublets(threshold=0.15)
            scrub.plot_histogram()[0].savefig(args.key + "_scrublet_hist.png")
            all_scores = numpy.vstack((doublet_scores, predicted_doublets))
            all_scores = numpy.transpose(all_scores)
            numpy.savetxt(args.key + "_scrublet_out.csv", all_scores, delimiter=",")
    else:
        filename = args.key + "_scrublet_hist.png"
        image = Image.new(mode = "RGB", size = (800,300), color = "white")
        draw = ImageDraw.Draw(image)
        draw.text((120,140), "Scrublet skipped by request.", fill = "black")
        image.save(filename)
        f = open(args.key + "_scrublet_out.csv", 'w')
        f.write("scrublet skipped")
        f.close()
