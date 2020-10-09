import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


class USGSFigure:
    def __init__(
            self, figure_type="map", family="Arial Narrow", font_path=None,
            verbose=False
    ):
        """Create a USGSFigure object

        Parameters
        ----------
        figure_type : str
            figure type ("map", "graph")
        family : str
            font family name (default is Arial Narrow)
        verbose : bool
            boolean that define if debug information should be written
        """
        # initialize members
        self.family = None
        self.figure_type = None
        self.verbose = verbose

        self.set_font_family(family=family, font_path=font_path)
        self.set_specifications(figure_type=figure_type)

    def set_specifications(self, figure_type="map"):
        """Set matplotlib parameters

        Parameters
        ----------
        figure_type : str
            figure type ("map", "graph")

        Returns
        -------

        """
        self.figure_type = self._validate_figure_type(figure_type)

    def set_font_family(self, family="Arial Narrow", font_path=None):
        """Set font family

        Parameters
        ----------
        family : str
            font family (default is Arial Narrow)
        font_path : str
            path to fonts not available to matplotlib (not implemented)

        Returns
        -------

        """
        if font_path is not None:
            errmsg = "specification of font_path is not implemented"
            raise NotImplemented(errmsg)
        self.family = self._set_fontfamily(family)

    def graph_legend(self, ax=None, handles=None, labels=None, **kwargs):
        """Add a USGS-style legend to a matplotlib axis object

        Parameters
        ----------
        ax : axis object
            matplotlib axis object (default is None)
        handles : list
            list of legend handles
        labels : list
            list of labels for legend handles
        kwargs : kwargs
            matplotlib legend kwargs

        Returns
        -------
        leg : object
            matplotlib legend object

        """
        if ax is None:
            ax = plt.gca()

        font = self._set_fontspec(bold=True, italic=False)
        if handles is None or labels is None:
            handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles, labels, prop=font, **kwargs)

        # add title to legend
        if "title" in kwargs:
            title = kwargs.pop("title")
        else:
            title = None
        leg = self.graph_legend_title(leg, title=title)
        return leg

    def graph_legend_title(self, leg, title=None):
        """Set the legend title for a matplotlib legend object

        Parameters
        ----------
        leg : legend object
            matplotlib legend object
        title : str
            title for legend

        Returns
        -------
        leg : object
            matplotlib legend object

        """
        if title is None:
            title = "EXPLANATION"
        elif title.lower() == "none":
            title = None
        font = self._set_fontspec(bold=True, italic=False)
        leg.set_title(title, prop=font)
        return leg

    def heading(self, ax=None, letter=None, heading=None, x=0.00, y=1.01, idx=None):
        """Add a USGS-style heading to a matplotlib axis object

        Parameters
        ----------
        ax : axis object
            matplotlib axis object (default is None)
        letter : str
            string that defines the subplot (A, B, C, etc.)
        heading : str
            text string
        x : float
            location of the heading in the x-direction in normalized plot dimensions
            ranging from 0 to 1 (default is 0.00)
        y : float
            location of the heading in the y-direction in normalized plot dimensions
            ranging from 0 to 1 (default is 1.01)
        idx : int
            index for programatically generating the heading letter when letter
            is None and idx is not None. idx = 0 will generate A (default is None)

        Returns
        -------
        text : object
            matplotlib text object

        """
        if ax is None:
            ax = plt.gca()

        if letter is None and idx is not None:
            letter = chr(ord("A") + idx)

        text = None
        if letter is not None:
            font = self._set_fontspec(bold=True, italic=True)
            if heading is None:
                letter = letter.replace(".", "")
            else:
                letter = letter.rstrip()
                if letter[-1] != ".":
                    letter += "."
                letter += " "
            ax.text(
                x,
                y,
                letter,
                va="bottom",
                ha="left",
                fontdict=font,
                transform=ax.transAxes,
            )
            bbox = ax.get_window_extent().transformed(
                plt.gcf().dpi_scale_trans.inverted()
            )
            width = bbox.width * 25.4  # inches to mm
            x += len(letter) * 1.0 / width
        if heading is not None:
            font = self._set_fontspec(bold=True, italic=False)
            text = ax.text(
                x,
                y,
                heading,
                va="bottom",
                ha="left",
                fontdict=font,
                transform=ax.transAxes,
            )
        return text

    def add_text(
            self,
            ax=None,
            text="",
            x=0.0,
            y=0.0,
            transform=True,
            bold=True,
            italic=True,
            fontsize=9,
            ha="left",
            va="bottom",
            **kwargs
    ):
        """Add USGS-style text to a axis object

        Parameters
        ----------
        ax : axis object
            matplotlib axis object (default is None)
        text : str
            text string
        x : float
            x-location of text string (default is 0.)
        y : float
            y-location of text string (default is 0.)
        transform : bool
            boolean that determines if a transformed (True) or data (False) coordinate
            system is used to define the (x, y) location of the text string
            (default is True)
        bold : bool
            boolean indicating if bold font (default is True)
        italic : bool
            boolean indicating if italic font (default is True)
        fontsize : int
            font size (default is 9 points)
        ha : str
            matplotlib horizontal alignment keyword (default is left)
        va : str
            matplotlib vertical alignment keyword (default is bottom)
        kwargs : dict
            dictionary with valid matplotlib text object keywords

        Returns
        -------
        text_obj : object
            matplotlib text object

        """
        if ax is None:
            ax = plt.gca()

        if transform:
            transform = ax.transAxes
        else:
            transform = ax.transData

        font = self._set_fontspec(bold=bold, italic=italic, fontsize=fontsize)

        text_obj = ax.text(
            x, y, text, va=va, ha=ha, fontdict=font, transform=transform, **kwargs
        )
        return text_obj

    def add_annotation(
            self,
            ax=None,
            text="",
            xy=None,
            xytext=None,
            bold=True,
            italic=True,
            fontsize=9,
            ha="left",
            va="bottom",
            **kwargs
    ):
        """Add an annotation to a axis object

        Parameters
        ----------
        ax : axis object
            matplotlib axis object (default is None)
        text : str
            text string
        xy : tuple
            tuple with the location of the annotation (default is None)
        xytext : tuple
            tuple with the location of the text
        bold : bool
            boolean indicating if bold font (default is True)
        italic : bool
            boolean indicating if italic font (default is True)
        fontsize : int
            font size (default is 9 points)
        ha : str
            matplotlib horizontal alignment keyword (default is left)
        va : str
            matplotlib vertical alignment keyword (default is bottom)
        kwargs : dict
            dictionary with valid matplotlib annotation object keywords

        Returns
        -------
        ann_obj : object
            matplotlib annotation object

        """
        if ax is None:
            ax = plt.gca()

        if xy is None:
            xy = (0.0, 0.0)

        if xytext is None:
            xytext = (0.0, 0.0)

        font = self._set_fontspec(bold=bold, italic=italic, fontsize=fontsize)

        # add font information to kwargs
        if kwargs is None:
            kwargs = font
        else:
            for key, value in font.items():
                kwargs[key] = value

        # create annotation
        ann_obj = ax.annotate(text, xy, xytext, va=va, ha=ha, **kwargs)

        return ann_obj

    def remove_edge_ticks(self, ax=None):
        """Remove unnecessary ticks on the edges of the plot

        Parameters
        ----------
        ax : axis object
            matplotlib axis object (default is None)

        Returns
        -------
        ax : axis object
            matplotlib axis object

        """
        if ax is None:
            ax = plt.gca()

        # update tick objects
        plt.draw()

        # get min and max value and ticks
        ymin, ymax = ax.get_ylim()

        # check for condition where y-axis values are reversed
        if ymax < ymin:
            y = ymin
            ymin = ymax
            ymax = y
        yticks = ax.get_yticks()

        if self.verbose:
            print("y-axis: ", ymin, ymax)
            print(yticks)

        # remove edge ticks on y-axis
        ticks = ax.yaxis.majorTicks
        for iloc in [0, -1]:
            if np.allclose(float(yticks[iloc]), ymin):
                ticks[iloc].tick1line.set_visible = False
                ticks[iloc].tick2line.set_visible = False
            if np.allclose(float(yticks[iloc]), ymax):
                ticks[iloc].tick1line.set_visible = False
                ticks[iloc].tick2line.set_visible = False

        # get min and max value and ticks
        xmin, xmax = ax.get_xlim()

        # check for condition where x-axis values are reversed
        if xmax < xmin:
            x = xmin
            xmin = xmax
            xmax = x

        xticks = ax.get_xticks()
        if self.verbose:
            print("x-axis: ", xmin, xmax)
            print(xticks)

        # remove edge ticks on y-axis
        ticks = ax.xaxis.majorTicks
        for iloc in [0, -1]:
            if np.allclose(float(xticks[iloc]), xmin):
                ticks[iloc].tick1line.set_visible = False
                ticks[iloc].tick2line.set_visible = False
            if np.allclose(float(xticks[iloc]), xmax):
                ticks[iloc].tick1line.set_visible = False
                ticks[iloc].tick2line.set_visible = False

        return ax

    # protected methods
    # protected method
    def _validate_figure_type(self, figure_type):
        """Set figure type after validation of specified figure type

        Parameters
        ----------
        figure_type : str
            figure type ("map", "graph")

        Returns
        -------
        figure_type : str
            validated figure_type

        """
        # validate figure type
        valid_types = ("map", "graph")
        if figure_type not in valid_types:
            errmsg = "invalid figure_type specified ({}) ".format(
                figure_type
            ) + "valid types are '{}'.".format(", ".join(valid_types))
            raise ValueError(errmsg)

        # set figure_type
        if figure_type == "map":
            self._set_map_specifications()
        elif figure_type == "graph":
            self._set_map_specifications()

        return figure_type

    # protected method
    def _set_graph_specifications(self):
        """Set matplotlib rcparams to USGS-style specifications for graphs

        Returns
        -------

        """
        rc_dict = {
            "font.family": self.family,
            "font.size": 7,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "axes.linewidth": 0.5,
            "xtick.labelsize": 8,
            "xtick.top": True,
            "xtick.bottom": True,
            "xtick.major.size": 7.2,
            "xtick.minor.size": 3.6,
            "xtick.major.width": 0.5,
            "xtick.minor.width": 0.5,
            "xtick.direction": "in",
            "ytick.labelsize": 8,
            "ytick.left": True,
            "ytick.right": True,
            "ytick.major.size": 7.2,
            "ytick.minor.size": 3.6,
            "ytick.major.width": 0.5,
            "ytick.minor.width": 0.5,
            "ytick.direction": "in",
            "pdf.fonttype": 42,
            "savefig.dpi": 300,
            "savefig.transparent": True,
            "legend.fontsize": 9,
            "legend.frameon": False,
            "legend.markerscale": 1.0,
        }
        mpl.rcParams.update(rc_dict)

    # protected method
    def _set_map_specifications(self):
        """Set matplotlib rcparams to USGS-style specifications for maps

        Returns
        -------

        """
        rc_dict = {
            "font.family": self.family,
            "font.size": 7,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "axes.linewidth": 0.5,
            "xtick.labelsize": 7,
            "xtick.top": True,
            "xtick.bottom": True,
            "xtick.major.size": 7.2,
            "xtick.minor.size": 3.6,
            "xtick.major.width": 0.5,
            "xtick.minor.width": 0.5,
            "xtick.direction": "in",
            "ytick.labelsize": 7,
            "ytick.left": True,
            "ytick.right": True,
            "ytick.major.size": 7.2,
            "ytick.minor.size": 3.6,
            "ytick.major.width": 0.5,
            "ytick.minor.width": 0.5,
            "ytick.direction": "in",
            "pdf.fonttype": 42,
            "savefig.dpi": 300,
            "savefig.transparent": True,
            "legend.fontsize": 9,
            "legend.frameon": False,
            "legend.markerscale": 1.0,
        }
        mpl.rcParams.update(rc_dict)

    # protected method
    def _set_fontspec(self, bold=True, italic=True, fontsize=9):
        """Create fontspec dictionary for matplotlib pyplot objects

        Parameters
        ----------
        bold : bool
            boolean indicating if font is bold (default is True)
        italic : bool
            boolean indicating if font is italic (default is True)
        fontsize : int
            font size (default is 9 point)


        Returns
        -------

        """
        if "Univers" in self.family:
            reset_family = True
        else:
            reset_family = False
            family = self.family

        if bold:
            weight = "bold"
            if reset_family:
                family = "Univers 67"
        else:
            weight = "normal"
            if reset_family:
                family = "Univers 57"

        if italic:
            if reset_family:
                family += " Condensed Oblique"
                style = "oblique"
            else:
                style = "italic"
        else:
            if reset_family:
                family += " Condensed"
            style = "normal"

        # define fontspec dictionary
        fontspec = {
            "family": family,
            "size": fontsize,
            "weight": weight,
            "style": style,
        }

        if self.verbose:
            sys.stdout.write("font specifications:\n ")
            for key, value in fontspec.items():
                sys.stdout.write("{}={} ".format(key, value))
            sys.stdout.write("\n")

        return fontspec

    def _set_fontfamily(self, family):
        """Set font family to Liberation Sans Narrow on linux if default Arial Narrow
        is being used

        Parameters
        ----------
        family : str
            font family name (default is Arial Narrow)

        Returns
        -------
        family : str
            font family name

        """
        if sys.platform.lower() in ("linux",):
            if family == "Arial Narrow":
                family = "Liberation Sans Narrow"
        return family
