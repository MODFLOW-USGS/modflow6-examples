import os

def build_table(caption, fpth, arr, headings=None, col_widths=None):
    if headings is None:
        headings = arr.dtype.names
    ncols = len(arr.dtype.names)
    if not fpth.endswith(".tex"):
        fpth += ".tex"
    label = "tab:{}".format(os.path.basename(fpth).replace(".tex", ""))

    line = get_header(caption, label, headings, col_widths=col_widths)

    for idx in range(arr.shape[0]):
        if idx % 2 != 0:
            line += "\t\t\\rowcolor{Gray}\n"
        line += "\t\t"
        for jdx, name in enumerate(arr.dtype.names):
            line += "{}".format(arr[name][idx])
            if jdx < ncols - 1:
                line += " & "
        line += " \\\\\n"

    # footer
    line += get_footer()

    f = open(fpth, "w")
    f.write(line)
    f.close()

    return

def get_header(caption, label, headings, col_widths=None, center=True, firsthead=False):
    ncol = len(headings)
    if col_widths is None:
        dx = 0.8 / float(ncol)
        col_widths = [dx for idx in range(ncol)]
    if center:
        align = "p"
    else:
        align = "p"

    header = "\\small\n"
    header += "\\begin{longtable}[!htbp]{\n"
    for col_width in col_widths:
        header += 38 * " " + \
                  "{}".format(align) + \
                  "{{{}\\linewidth-2\\arraycolsep}}\n".format(col_width)
    header += 38 * " " + "}\n"
    header += "\t\\caption{{{}}} \\label{{{}}} \\\\\n\n".format(caption, label)

    if firsthead:
        header += "\t\\hline \\hline\n"
        header +=  "\t\\rowcolor{Gray}\n"
        header += "\t"
        for idx, s in enumerate(headings):
            header += "\\textbf{{{}}}".format(s)
            if idx < len(headings) - 1:
                header += " & "
        header += "  \\\\\n"
        header +=  "\t\\hline\n"
        header +=  "\t\\endfirsthead\n\n"

    header += "\t\\hline \\hline\n"
    header +=  "\t\\rowcolor{Gray}\n"
    header += "\t"
    for idx, s in enumerate(headings):
        header += "\\textbf{{{}}}".format(s)
        if idx < len(headings) - 1:
            header += " & "
    header += "  \\\\\n"
    header +=  "\t\\hline\n"
    header +=  "\t\\endhead\n\n"

    return header

def get_footer():
    return "\t\\hline \\hline\n\\end{longtable}\n\\normalsize\n\n"

def exp_format(v):
    s = "{:.2e}".format(v)
    s = s.replace("e-0", "e-")
    s = s.replace("e+0", "e+")
    # s = s.replace("e", " \\times 10^{") + "}$"
    return s

def float_format(v, fmt="{:.2f}"):
    return fmt.format(v)

def int_format(v):
    return "{:d}".format(v)