#!/usr/bin/env python3
"""
Grayscale print copies of the figure TIFFs (color stays for the online version).

For every color TIFF in a folder this writes a new file next to it, <name>_print_gray.tif[f]:
8-bit single-channel grayscale (ITU-R 601 luma, the conversion most journals and printers use), LZW-compressed,
with the same pixel size and the same dpi tag as the original (600). The color originals are never touched,
and files that already end in _print_gray are skipped, so it is safe to run again after re-rendering.

    python3 make_print_gray.py  path/to/figure            # every .tif / .tiff in the folder
    python3 make_print_gray.py  path/to/figure  Fig3 Fig4  # only files whose name contains one of these words

Requires: Pillow.
"""
import os
import sys
from PIL import Image

Image.MAX_IMAGE_PIXELS = None            # the Fig. 2 canvas is large (8643 x 5919 px)
SUFFIX = "_print_gray"


def gray_copy(path):
    base, ext = os.path.splitext(path)
    out = base + SUFFIX + ext
    im = Image.open(path)
    dpi = im.info.get("dpi", (600, 600))
    im.convert("L").save(out, compression="tiff_lzw", dpi=dpi)     # "L" = 8-bit grayscale (luma = 0.299 R + 0.587 G + 0.114 B)
    return out, im.size, dpi


if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    words = sys.argv[2:]
    names = sorted(f for f in os.listdir(folder)
                   if f.lower().endswith((".tif", ".tiff")) and SUFFIX not in f and (not words or any(w in f for w in words)))
    if not names:
        sys.exit(f"No color TIFFs found in {folder}")
    for f in names:
        out, size, dpi = gray_copy(os.path.join(folder, f))
        print(f"{f}  ->  {os.path.basename(out)}   {size[0]} x {size[1]} px, dpi {tuple(round(float(d)) for d in dpi)}")
