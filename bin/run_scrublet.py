#!/usr/bin/env python
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
    parser.add_argument('--mat', required=True, help='input matrix.')
    args = parser.parse_args()

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
        image = Image.new(mode = "RGB", size = (250,50), color = "white")
        draw = ImageDraw.Draw(image)
        draw.text((10,10), "Scrublet failed. This is generally \nbecause there aren't enough cells.", fill = "black")
        image.save(filename)
        numpy.savetxt(args.key + "_scrublet_out.csv", all_scores, fmt="%s", delimiter=",")
    except (AttributeError):
        predicted_doublets = scrub.call_doublets(threshold=0.15)
        scrub.plot_histogram()[0].savefig(args.key + "_scrublet_hist.png")
        all_scores = numpy.vstack((doublet_scores, predicted_doublets))
        all_scores = numpy.transpose(all_scores)
        numpy.savetxt(args.key + "_scrublet_out.csv", all_scores, delimiter=",")
