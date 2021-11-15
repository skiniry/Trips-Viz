import mpld3
from mpld3 import plugins,utils


import collections
import json
import uuid
import matplotlib








class PointLabelTooltip(plugins.PluginBase):
    """A Plugin to enable a tooltip: text which hovers over points.
    Parameters
    ----------
    points : matplotlib Collection or Line2D object
        The figure element to apply the tooltip to
    labels : array or None
        If supplied, specify the labels for each point in points.  If not
        supplied, the (x, y) values will be used.
    hoffset, voffset : integer
        The number of pixels to offset the tooltip text.  Default is
        hoffset = 0, voffset = 10
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from mpld3 import fig_to_html, plugins
    >>> fig, ax = plt.subplots()
    >>> points = ax.plot(range(10), 'o')
    >>> plugins.connect(fig, PointLabelTooltip(points[0]))
    >>> fig_to_html(fig)
    """
    def __init__(self, points, labels=None,
                 hoffset=0, voffset=10, location="mouse"):
        if location not in ["bottom left", "top left", "bottom right",
                            "top right", "mouse"]:
            raise ValueError("invalid location: {0}".format(location))
        if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None
        self.dict_ = {"type": "tooltip",
                      "id": plugins.get_id(points, suffix),
                      "labels": labels,
                      "hoffset": hoffset,
                      "voffset": voffset,
                      "location": location}


class LineLabelTooltip(plugins.PluginBase):
    """A Plugin to enable a tooltip: text which hovers over a line.
    Parameters
    ----------
    line : matplotlib Line2D object
        The figure element to apply the tooltip to
    label : string
        If supplied, specify the labels for each point in points.  If not
        supplied, the (x, y) values will be used.
    hoffset, voffset : integer
        The number of pixels to offset the tooltip text.  Default is
        hoffset = 0, voffset = 10
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from mpld3 import fig_to_html, plugins
    >>> fig, ax = plt.subplots()
    >>> lines = ax.plot(range(10), 'o')
    >>> plugins.connect(fig, LineLabelTooltip(lines[0]))
    >>> fig_to_html(fig)
    """
    def __init__(self, points, label=None,
                 hoffset=0, voffset=10, location="mouse"):
        if location not in ["bottom left", "top left", "bottom right",
                            "top right", "mouse"]:
            raise ValueError("invalid location: {0}".format(location))
        self.dict_ = {"type": "tooltip",
                      "id": plugins.get_id(points),
                      "labels": label if label is None else [label],
                      "hoffset": hoffset,
                      "voffset": voffset,
                      "location": location}


class LinkedBrush(plugins.PluginBase):
    """A Plugin to enable linked brushing between plots
    Parameters
    ----------
    points : matplotlib Collection or Line2D object
        A representative of the scatter plot elements to brush.
    button : boolean, optional
        if True (default), then add a button to enable/disable zoom behavior
    enabled : boolean, optional
        specify whether the zoom should be enabled by default. default=True.
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from mpld3 import fig_to_html, plugins
    >>> X = np.random.random((3, 100))
    >>> fig, ax = plt.subplots(3, 3)
    >>> for i in range(2):
    ...     for j in range(2):
    ...         points = ax[i, j].scatter(X[i], X[j])
    >>> plugins.connect(fig, LinkedBrush(points))
    >>> fig_to_html(fig)
    Notes
    -----
    Notice that in the above example, only one of the four sets of points is
    passed to the plugin. This is all that is needed: for the sake of efficient
    data storage, mpld3 keeps track of which plot objects draw from the same
    data.
    Also note that for the linked brushing to work correctly, the data must
    not contain any NaNs. The presence of NaNs makes the different data views
    have different sizes, so that mpld3 is unable to link the related points.
    """

    def __init__(self, points, button=True, enabled=True):
        if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None

        self.dict_ = {"type": "linkedbrush",
                      "button": button,
                      "enabled": enabled,
                      "id": plugins.get_id(points, suffix)}


class PointHTMLTooltip(plugins.PluginBase):
    """A Plugin to enable an HTML tooltip:
    formated text which hovers over points.
    Parameters
    ----------
    points : matplotlib Collection or Line2D object
        The figure element to apply the tooltip to
    labels : list
        The labels for each point in points, as strings of unescaped HTML.
    hoffset, voffset : integer, optional
        The number of pixels to offset the tooltip text.  Default is
        hoffset = 0, voffset = 10
    css : str, optional
        css to be included, for styling the label html if desired
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from mpld3 import fig_to_html, plugins
    >>> fig, ax = plt.subplots()
    >>> points = ax.plot(range(10), 'o')
    >>> labels = ['<h1>{title}</h1>'.format(title=i) for i in range(10)]
    >>> plugins.connect(fig, PointHTMLTooltip(points[0], labels))
    >>> fig_to_html(fig)
    """

    JAVASCRIPT = """
    mpld3.register_plugin("htmltooltip", HtmlTooltipPlugin);
    HtmlTooltipPlugin.prototype = Object.create(mpld3.Plugin.prototype);
    HtmlTooltipPlugin.prototype.constructor = HtmlTooltipPlugin;
    HtmlTooltipPlugin.prototype.requiredProps = ["id"];
    HtmlTooltipPlugin.prototype.defaultProps = {labels:null,
                                                hoffset:0,
                                                voffset:10};
    function HtmlTooltipPlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    HtmlTooltipPlugin.prototype.draw = function(){
       var obj = mpld3.get_element(this.props.id);
       var labels = this.props.labels;
       var tooltip = d3.select("body").append("div")
                    .attr("class", "mpld3-tooltip")
                    .style("position", "absolute")
                    .style("z-index", "10")
                    .style("visibility", "hidden");
       obj.elements()
           .on("mouseover", function(d, i){
                              tooltip.html(labels[i])
                                     .style("visibility", "visible");})
           .on("mousemove", function(d, i){
                  tooltip
                    .style("top", d3.event.pageY + this.props.voffset + "px")
                    .style("left",d3.event.pageX + this.props.hoffset + "px");
                 }.bind(this))
           .on("mouseout",  function(d, i){
                           tooltip.style("visibility", "hidden");});
    };
    """

    def __init__(self, points, labels=None,
                 hoffset=0, voffset=10, css=None):
        self.points = points
        self.labels = labels
        self.voffset = voffset
        self.hoffset = hoffset
        self.css_ = css or ""
        if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None
        self.dict_ = {"type": "htmltooltip",
                      "id": plugins.get_id(points, suffix),
                      "labels": labels,
                      "hoffset": hoffset,
                      "voffset": voffset}


class LineHTMLTooltip(plugins.PluginBase):
    """A Plugin to enable an HTML tooltip:
    formated text which hovers over points.
    Parameters
    ----------
    points : matplotlib Line2D object
        The figure element to apply the tooltip to
    label : string
        The label for the line, as strings of unescaped HTML.
    hoffset, voffset : integer, optional
        The number of pixels to offset the tooltip text.  Default is
        hoffset = 0, voffset = 10
    css : str, optional
        css to be included, for styling the label html if desired
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from mpld3 import fig_to_html, plugins
    >>> fig, ax = plt.subplots()
    >>> lines = ax.plot(range(10))
    >>> label = '<h1>line {title}</h1>'.format(title='A')
    >>> plugins.connect(fig, LineHTMLTooltip(lines[0], label))
    >>> fig_to_html(fig)
    """

    JAVASCRIPT = """
    mpld3.register_plugin("linehtmltooltip", LineHTMLTooltip);
    LineHTMLTooltip.prototype = Object.create(mpld3.Plugin.prototype);
    LineHTMLTooltip.prototype.constructor = LineHTMLTooltip;
    LineHTMLTooltip.prototype.requiredProps = ["id"];
    LineHTMLTooltip.prototype.defaultProps = {label:null,
                                              hoffset:0,
                                              voffset:10};
    function LineHTMLTooltip(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    LineHTMLTooltip.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id, this.fig);
        var label = this.props.label
        var tooltip = d3.select("body").append("div")
                    .attr("class", "mpld3-tooltip")
                    .style("position", "absolute")
                    .style("z-index", "10")
                    .style("visibility", "hidden");
        obj.elements()
           .on("mouseover", function(d, i){
                               tooltip.html(label)
                                      .style("visibility", "visible");
                                     })
            .on("mousemove", function(d, i){
                  tooltip
                    .style("top", d3.event.pageY + this.props.voffset + "px")
                    .style("left",d3.event.pageX + this.props.hoffset + "px");
                 }.bind(this))
           .on("mouseout",  function(d, i){
                           tooltip.style("visibility", "hidden");})
    };
    """

    def __init__(self, line, label=None,
                 hoffset=0, voffset=10,
                 css=None):
        self.line = line
        self.label = label
        self.voffset = voffset
        self.hoffset = hoffset
        self.css_ = css or ""
        self.dict_ = {"type": "linehtmltooltip",
                      "id": plugins.get_id(line),
                      "label": label,
                      "hoffset": hoffset,
                      "voffset": voffset}


class InteractiveLegendPlugin(plugins.PluginBase):
	"""A plugin for an interactive legends.

	Inspired by http://bl.ocks.org/simzou/6439398

	Parameters
	----------
	plot_elements : iterable of matplotliblib elements
		the elements to associate with a given legend items
	labels : iterable of strings
		The labels for each legend element
	ax :  matplotlib axes instance, optional
		the ax to which the legend belongs. Default is the first
		axes. The legend will be plotted to the right of the specified
		axes
	alpha_sel : float, optional
		the alpha value to apply to the plot_element(s) associated
		with the legend item when the legend item is selected.
		Default is 1.0
	alpha_unsel : float, optional
		the alpha value to apply to the plot_element(s) associated
		with the legend item when the legend item is unselected.
		Default is 0.2
	xoffset : int, optional
		apply x offset to the rectangles and labels
	yoffset : int, optional
		apply y offset to the rectangles and labels
	start_visible : boolean, optional (could be a list of booleans)
		defines if objects should start selected on not.
	fontsize: int, optional
		defines font size of the legend labels
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld4 import fig_to_html, plugins
	>>> N_paths = 5
	>>> N_steps = 100
	>>> x = np.linspace(0, 10, 100)
	>>> y = 0.1 * (np.random.random((N_paths, N_steps)) - 0.5)
	>>> y = y.cumsum(1)
	>>> fig, ax = plt.subplots()
	>>> labels = ["a", "b", "c", "d", "e"]
	>>> line_collections = ax.plot(x, y.T, lw=4, alpha=0.1)
	>>> interactive_legend = plugins.InteractiveLegendPlugin(line_collections,
	...                                                      labels,
	...                                                      alpha_unsel=0.1)
	>>> plugins.connect(fig, interactive_legend)
	>>> fig_to_html(fig)
	"""

	JAVASCRIPT = """
	mpld3.register_plugin("interactive_legend", InteractiveLegend);
	InteractiveLegend.prototype = Object.create(mpld3.Plugin.prototype);
	InteractiveLegend.prototype.constructor = InteractiveLegend;
	InteractiveLegend.prototype.requiredProps = ["element_ids", "labels"];
	InteractiveLegend.prototype.defaultProps = {"ax":null,
												"alpha_sel":1.0,
												"alpha_unsel":0,
												"xoffset":0,
												"yoffset":0,
												"fontsize":17,
												"start_visible":[false,false,false,false,false,false,false,false,false,false,false]}
	function InteractiveLegend(fig, props){
		mpld3.Plugin.call(this, fig, props);
	};

	InteractiveLegend.prototype.draw = function(){
		console.log(this);
		var alpha_sel = this.props.alpha_sel;
		var alpha_unsel = this.props.alpha_unsel;
		var xoffset = this.props.xoffset;
		var yoffset = this.props.yoffset;
		var fontsize = this.props.fontsize;

		var legendItems = new Array();
		for(var i=0; i<this.props.labels.length; i++){
			var obj = {};
			obj.label = this.props.labels[i];

			var element_id = this.props.element_ids[i];
			mpld3_elements = [];
			for(var j=0; j<element_id.length; j++){
				var mpld3_element = mpld3.get_element(element_id[j], this.fig);

				// mpld3_element might be null in case of Line2D instances
				// for we pass the id for both the line and the markers. Either
				// one might not exist on the D3 side
				if(mpld3_element){
					mpld3_element.visible = this.props.start_visible[i];
					mpld3_elements.push(mpld3_element);
				}
			}

			obj.mpld3_elements = mpld3_elements;
			
			obj.visible = this.props.start_visible[i]; // should become be setable from python side
			legendItems.push(obj);

		}
		console.log("LEGEND ITEMS:"+legendItems);

		// determine the axes with which this legend is associated
		var ax = this.props.ax
		if(!ax){
			ax = this.fig.axes[0];
		} else{
			ax = mpld3.get_element(ax, this.fig);
		}

		// add a legend group to the canvas of the figure
		var legend = this.fig.canvas.append("svg:g")
							   .attr("class", "legend");

		// add the cds_areagles
		legend.selectAll("rect")
				.data(legendItems)
			 .enter().append("rect")
				.attr("height",10)
				.attr("width", 20)
				.attr("x",ax.width+1+ax.position[0]-(xoffset))
				.attr("y",function(d,i) {
							return ax.position[1]+ i * 25 + (yoffset);}) // should be -10
				.attr("stroke", get_color)
				.attr("class", "legend-box")
				.style("fill", function(d, i) {
							return d.visible ? get_color(d) : "white";})
				.on("click", click);

		// add the labels
		legend.selectAll("text")
			  .data(legendItems)
			  .enter().append("text")
			  .attr("x", function (d) {
							return ax.width+5+ax.position[0] - (xoffset-25);}) // should be +25
			  .attr("y", function(d,i) {
							return ax.position[1]+ (i * 25)+(yoffset+10);})
			  .attr("font-size",fontsize+"px")
			  .text(function(d) { return d.label });

		// specify the action on click
		function click(d,i){
			d.visible = !d.visible;
			d3.select(this)
			  .style("fill",function(d, i) {
				return d.visible ? get_color(d) : "white";
			  })

			for(var i=0; i<d.mpld3_elements.length; i++){
				var type = d.mpld3_elements[i].constructor.name;
				if(type =="mpld3_Line"){
					d3.select(d.mpld3_elements[i].path[0][0])
						.style("stroke-opacity",
								d.visible ? alpha_sel : alpha_unsel);
				} else if((type=="mpld3_PathCollection")||
						 (type=="mpld3_Markers")){
					d3.selectAll(d.mpld3_elements[i].pathsobj[0])
						.style("stroke-opacity",
								d.visible ? alpha_sel : alpha_unsel)
						.style("fill-opacity",
								d.visible ? alpha_sel : alpha_unsel);
				} else if(type=="mpld3_Path"){
					var current_alpha = d.mpld3_elements[i].props.alpha;
					var current_alpha_unsel = current_alpha * alpha_unsel;
					var current_alpha_over = current_alpha * alpha_sel;
					d3.select(d.mpld3_elements[i].path[0][0])
						.style("stroke-opacity",
												d.visible ? alpha_sel : current_alpha_unsel)
						.style("fill-opacity",
												d.visible ? alpha_sel : current_alpha_unsel)
				}else{
					console.log(type + " not yet supported");
				}
			}
		};

		// helper function for determining the color of the cds_areagles
		function get_color(d){
			var type = d.mpld3_elements[0].constructor.name;
			var color = "black";
			if(type =="mpld3_Line"){
				color = d.mpld3_elements[0].props.edgecolor;
			} else if((type=="mpld3_PathCollection")||
					 (type=="mpld3_Markers")){
				color = d.mpld3_elements[0].props.facecolors[0];
			} else if(type=="mpld3_Path"){
				color = d.mpld3_elements[0].props.facecolor;
			}else{
				console.log(type + " not yet supported");
			}
			return color;
		};
	};
	"""

	css_ = """
	.legend-box {
	  cursor: pointer;
	}
	"""

	def __init__(self, plot_elements, labels, ax=None,
				 alpha_sel=1, alpha_unsel=0.2, xoffset=0, yoffset=0,start_visible=False,fontsize=18):

		self.ax = ax

		if ax:
			ax = plugins.get_id(ax)

		# start_visible could be a list
		if isinstance(start_visible, bool):
			start_visible = [start_visible] * len(labels)
		elif not len(start_visible) == len(labels):
			raise ValueError("{} out of {} visible params has been set".format(len(start_visible), len(labels)))

		mpld4_element_ids = self._determine_mpld4ids(plot_elements)
		self.mpld4_element_ids = mpld4_element_ids
		self.dict_ = {"type": "interactive_legend",
					  "element_ids": mpld4_element_ids,
					  "labels": labels,
					  "ax": ax,
					  "alpha_sel": alpha_sel,
					  "alpha_unsel": alpha_unsel,
					  "xoffset":xoffset,
					  "yoffset":yoffset,
					  "start_visible":start_visible,
					  "fontsize":fontsize}

	def _determine_mpld4ids(self, plot_elements):
		"""
		Helper function to get the mpld4_id for each
		of the specified elements.
		"""
		mpld4_element_ids = []

		# There are two things being done here. First,
		# we make sure that we have a list of lists, where
		# each inner list is associated with a single legend
		# item. Second, in case of Line2D object we pass
		# the id for both the marker and the line.
		# on the javascript side we filter out the nulls in
		# case either the line or the marker has no equivalent
		# D3 representation.
		for entry in plot_elements:
			ids = []
			if isinstance(entry, collections.Iterable):
				for element in entry:
					mpld4_id = plugins.get_id(element)
					ids.append(mpld4_id)
					if isinstance(element, matplotlib.lines.Line2D):
						mpld4_id = plugins.get_id(element, 'pts')
						ids.append(mpld4_id)
			else:
				ids.append(plugins.get_id(entry))
				if isinstance(entry, matplotlib.lines.Line2D):
					mpld4_id = plugins.get_id(entry, 'pts')
					ids.append(mpld4_id)
			mpld4_element_ids.append(ids)
		return mpld4_element_ids


class PointClickableHTMLTooltip(plugins.PluginBase):
    """A plugin for pop-up windows with data with rich HTML
    Parameters
    ----------
    points : matplotlib Collection object
        The figure element to apply the tooltip to
    labels : list
        The labels for each point in points, as strings of unescaped HTML.
    targets : list
        The target data or rich HTML to be displayed when each collection element is clicked
    hoffset, voffset : integer, optional
        The number of pixels to offset the tooltip text.  Default is
        hoffset = 0, voffset = 10
    css : str, optional
        css to be included, for styling the label html and target data/tables, if desired
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from mpld3 import plugins
    >>> fig, ax = plt.subplots(1,1)
    >>> xx = yy = range(10)
    >>> scat = ax.scatter(xx, range(10))
    >>> targets = map(lambda (x, y): "<marquee>It works!<br><h1>{}, {}</h1></marquee>".format(x, y),
    >>>               zip(xx, yy))
    >>> labels = map(lambda (x, y): "{}, {}".format(x,y), zip(xx, yy))
    >>> from mpld3.plugins import PointClickableHTMLTooltip
    >>> plugins.connect(fig, PointClickableHTMLTooltip(scat, labels=labels, targets=targets))
    """

    JAVASCRIPT = """
    mpld3.register_plugin("clickablehtmltooltip", PointClickableHTMLTooltip);
    PointClickableHTMLTooltip.prototype = Object.create(mpld3.Plugin.prototype);
    PointClickableHTMLTooltip.prototype.constructor = PointClickableHTMLTooltip;
    PointClickableHTMLTooltip.prototype.requiredProps = ["id"];
    PointClickableHTMLTooltip.prototype.defaultProps = {labels:null,
                                                 targets:null,
                                                 hoffset:0,
                                                 voffset:10};
    function PointClickableHTMLTooltip(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    PointClickableHTMLTooltip.prototype.draw = function(){
       var obj = mpld3.get_element(this.props.id);
       var labels = this.props.labels;
       var targets = this.props.targets;
       var tooltip = d3.select("body").append("div")
                    .attr("class", "mpld3-tooltip")
                    .style("position", "absolute")
                    .style("z-index", "10")
                    .style("visibility", "hidden");
       obj.elements()
           .on("mouseover", function(d, i){
                  if ($(obj.elements()[0][0]).css( "fill-opacity" ) > 0 || $(obj.elements()[0][0]).css( "stroke-opacity" ) > 0) {
                              tooltip.html(labels[i])
                                     .style("visibility", "visible");
                              } })
           .on("mousedown", function(d, i){
                              window.open().document.write(targets[i]);
                               })
           .on("mousemove", function(d, i){
                  tooltip
                    .style("top", d3.event.pageY + this.props.voffset + "px")
                    .style("left",d3.event.pageX + this.props.hoffset + "px");
                 }.bind(this))
           .on("mouseout",  function(d, i){
                           tooltip.style("visibility", "hidden");});
    };
    """

    def __init__(self, points, labels=None, targets=None,
                 hoffset=2, voffset=-6, css=None):
        self.points = points
        self.labels = labels
        self.targets = targets
        self.voffset = voffset
        self.hoffset = hoffset
        self.css_ = css or ""
        if targets is not None:
            styled_targets = map(lambda x: self.css_ + x, targets)
        else:
            styled_targets = None


        if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None
        self.dict_ = {"type": "clickablehtmltooltip",
                      "id": plugins.get_id(points, suffix),
                      "labels": labels,
                      "targets": styled_targets,
                      "hoffset": hoffset,
                      "voffset": voffset}


class MouseXPosition(plugins.PluginBase):
    """Like MousePosition, but only show the X coordinate"""

    JAVASCRIPT = """
  mpld3.register_plugin("mousexposition", MouseXPositionPlugin);
  MouseXPositionPlugin.prototype = Object.create(mpld3.Plugin.prototype);
  MouseXPositionPlugin.prototype.constructor = MouseXPositionPlugin;
  MouseXPositionPlugin.prototype.requiredProps = [];
  MouseXPositionPlugin.prototype.defaultProps = {
    fontsize: 12,
    fmt: "0d"
  };
  function MouseXPositionPlugin(fig, props) {
    mpld3.Plugin.call(this, fig, props);
  }
  MouseXPositionPlugin.prototype.draw = function() {
    var fig = this.fig;
    var fmt = d3.format(this.props.fmt);
    var coords = fig.canvas.append("text").attr("class", "mpld3-coordinates").style("text-anchor", "end").style("font-size", this.props.fontsize).attr("x", this.fig.width - 5).attr("y", this.fig.height - 5);
    for (var i = 0; i < this.fig.axes.length; i++) {
      var update_coords = function() {
        var ax = fig.axes[i];
        return function() {
          var pos = d3.mouse(this), x = ax.x.invert(pos[0]), y = ax.y.invert(pos[1]);
          coords.text(fmt(x));
        };
      }();
      fig.axes[i].baseaxes.on("mousemove", update_coords).on("mouseout", function() {
        coords.text("");
      });
    }
  };"""
    """A Plugin to display coordinates for the current mouse position
    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>> from mpld3 import fig_to_html, plugins
    >>> fig, ax = plt.subplots()
    >>> points = ax.plot(range(10), 'o')
    >>> plugins.connect(fig, plugins.MouseXPosition())
    >>> fig_to_html(fig)
    """
    def __init__(self, fontsize=12, fmt="8.0f"):
        self.dict_ = {"type": "mousexposition",
                      "fontsize": fontsize,
                      "fmt": fmt}
        
        
class TopToolbar(plugins.PluginBase):
	"""Plugin for moving toolbar to top of figure"""

	JAVASCRIPT = """
	mpld3.register_plugin("toptoolbar", TopToolbar);
	TopToolbar.prototype = Object.create(mpld3.Plugin.prototype);
	TopToolbar.prototype.constructor = TopToolbar;
	TopToolbar.prototype.defaultProps = {"xoffset":0,
										 "yoffset":0}
	function TopToolbar(fig, props){
		mpld3.Plugin.call(this, fig, props);
	};

	TopToolbar.prototype.draw = function(){
	  // the toolbar svg doesn't exist
	  // yet, so first draw it
	  this.fig.toolbar.draw();
	  var xoffset = this.props.xoffset;
	  var yoffset = this.props.yoffset;

	  // set the x and y position of the toolbar
	  this.fig.toolbar.toolbar.attr("y", 800+yoffset);
	  this.fig.toolbar.toolbar.attr("x", 910+xoffset);

	  // set the width and height of the toolbar
	  this.fig.toolbar.toolbar.attr("width", 280);
	  this.fig.toolbar.toolbar.attr("height", 30);

	  // set the width and height of the buttons
	  this.fig.toolbar.buttonsobj.attr("width",24);
	  this.fig.toolbar.buttonsobj.attr("height",24);

	  // set the space between the buttons
	  this.fig.toolbar.buttonsobj.attr("x", function(d, i) {return i * 30;});

	  // show toolbar
	  this.fig.toolbar.buttonsobj.transition(0).attr("y", 0);

	  // remove event triggers so that toolbar doesnt dissapear when mouse leaves plot.
	  this.fig.canvas
		.on("mouseenter", null)
		.on("mouseleave", null)
		.on("touchenter", null)
		.on("touchstart", null);

	  // then remove the draw function,
	  // so that it is not called again
	  this.fig.toolbar.draw = function() {}
	}
	"""
	def __init__(self,xoffset=0, yoffset=0):
		self.dict_ = {"type": "toptoolbar",
					  "xoffset":xoffset,
					  "yoffset":yoffset} 
        

class DownloadProfile(plugins.PluginBase):
	"""A Plugin to enable an HTML tooltip:
	formated text which hovers over points.

	Parameters
	----------
	points : matplotlib Collection or Line2D object
		The figure element to apply the tooltip to
	labels : list
		The labels for each point in points, as strings of unescaped HTML.
	hoffset, voffset : integer, optional
		The number of pixels to offset the tooltip text.  Default is
		hoffset = 0, voffset = 10
	css : str, optional
		css to be included, for styling the label html if desired
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> points = ax.plot(range(10), 'o')
	>>> labels = ['<h1>{title}</h1>'.format(title=i) for i in range(10)]
	>>> plugins.connect(fig, PointHTMLTooltip(points[0], labels))
	>>> fig_to_html(fig)
	"""

	JAVASCRIPT = """




	/* FileSaver.js
	 * A saveAs() FileSaver implementation.
	 * 1.3.8
	 * 2018-03-22 14:03:47
	 *
	 * By Eli Grey, https://eligrey.com
	 * License: MIT
	 *   See https://github.com/eligrey/FileSaver.js/blob/master/LICENSE.md
	 */

	/*global self */
	/*jslint bitwise: true, indent: 4, laxbreak: true, laxcomma: true, smarttabs: true, plusplus: true */

	/*! @source http://purl.eligrey.com/github/FileSaver.js/blob/master/src/FileSaver.js */

	var saveAs = saveAs || (function(view) {
	"use strict";
	// IE <10 is explicitly unsupported
	if (typeof view === "undefined" || typeof navigator !== "undefined" && /MSIE [1-9]\./.test(navigator.userAgent)) {
		return;
	}
	var
	doc = view.document
	// only get URL when necessary in case Blob.js hasn't overridden it yet
	, get_URL = function() {
		return view.URL || view.webkitURL || view;
	}
	, save_link = doc.createElementNS("http://www.w3.org/1999/xhtml", "a")
	, can_use_save_link = "download" in save_link
	, click = function(node) {
		var event = new MouseEvent("click");
		node.dispatchEvent(event);
	}
	, is_safari = /constructor/i.test(view.HTMLElement) || view.safari
	, is_chrome_ios =/CriOS\/[\d]+/.test(navigator.userAgent)
	, setImmediate = view.setImmediate || view.setTimeout
	, throw_outside = function(ex) {
		setImmediate(function() {
			throw ex;
		}, 0);
	}
	, force_saveable_type = "application/octet-stream"
	// the Blob API is fundamentally broken as there is no "downloadfinished" event to subscribe to
	, arbitrary_revoke_timeout = 1000 * 40 // in ms
	, revoke = function(file) {
		var revoker = function() {
			if (typeof file === "string") { // file is an object URL
				get_URL().revokeObjectURL(file);
			} else { // file is a File
				file.remove();
			}
		};
		setTimeout(revoker, arbitrary_revoke_timeout);
	}
	, dispatch = function(filesaver, event_types, event) {
		event_types = [].concat(event_types);
		var i = event_types.length;
		while (i--) {
			var listener = filesaver["on" + event_types[i]];
			if (typeof listener === "function") {
				try {
					listener.call(filesaver, event || filesaver);
				} catch (ex) {
					throw_outside(ex);
				}
			}
		}
	}
	, auto_bom = function(blob) {
		// prepend BOM for UTF-8 XML and text/* types (including HTML)
		// note: your browser will automatically convert UTF-16 U+FEFF to EF BB BF
		if (/^\s*(?:text\/\S*|application\/xml|\S*\/\S*\+xml)\s*;.*charset\s*=\s*utf-8/i.test(blob.type)) {
			return new Blob([String.fromCharCode(0xFEFF), blob], {type: blob.type});
		}
		return blob;
	}
	, FileSaver = function(blob, name, no_auto_bom) {
		if (!no_auto_bom) {
			blob = auto_bom(blob);
		}
		// First try a.download, then web filesystem, then object URLs
		var
		filesaver = this
		, type = blob.type
		, force = type === force_saveable_type
		, object_url
		, dispatch_all = function() {
			dispatch(filesaver, "writestart progress write writeend".split(" "));
		}
		// on any filesys errors revert to saving with object URLs
		, fs_error = function() {
			if ((is_chrome_ios || (force && is_safari)) && view.FileReader) {
				// Safari doesn't allow downloading of blob urls
				var reader = new FileReader();
				reader.onloadend = function() {
					var url = is_chrome_ios ? reader.result : reader.result.replace(/^data:[^;]*;/, 'data:attachment/file;');
					var popup = view.open(url, '_blank');
					if(!popup) view.location.href = url;
							   url=undefined; // release reference before dispatching
							   filesaver.readyState = filesaver.DONE;
					dispatch_all();
				};
				reader.readAsDataURL(blob);
				filesaver.readyState = filesaver.INIT;
				return;
			}
			// don't create more object URLs than needed
			if (!object_url) {
				object_url = get_URL().createObjectURL(blob);
			}
			if (force) {
				view.location.href = object_url;
			} else {
				var opened = view.open(object_url, "_blank");
				if (!opened) {
					// Apple does not allow window.open, see https://developer.apple.com/library/safari/documentation/Tools/Conceptual/SafariExtensionGuide/WorkingwithWindowsandTabs/WorkingwithWindowsandTabs.html
					view.location.href = object_url;
				}
			}
			filesaver.readyState = filesaver.DONE;
			dispatch_all();
			revoke(object_url);
		}
		;
		filesaver.readyState = filesaver.INIT;

		if (can_use_save_link) {
			object_url = get_URL().createObjectURL(blob);
			setImmediate(function() {
				save_link.href = object_url;
				save_link.download = name;
				click(save_link);
				dispatch_all();
				revoke(object_url);
				filesaver.readyState = filesaver.DONE;
			}, 0);
			return;
		}

		fs_error();
	}
	, FS_proto = FileSaver.prototype
	, saveAs = function(blob, name, no_auto_bom) {
		return new FileSaver(blob, name || blob.name || "download", no_auto_bom);
	}
	;

	// IE 10+ (native saveAs)
	if (typeof navigator !== "undefined" && navigator.msSaveOrOpenBlob) {
		return function(blob, name, no_auto_bom) {
			name = name || blob.name || "download";

			if (!no_auto_bom) {
				blob = auto_bom(blob);
			}
			return navigator.msSaveOrOpenBlob(blob, name);
		};
	}

	// todo: detect chrome extensions & packaged apps
	//save_link.target = "_blank";

	FS_proto.abort = function(){};
	FS_proto.readyState = FS_proto.INIT = 0;
	FS_proto.WRITING = 1;
	FS_proto.DONE = 2;

	FS_proto.error =
	FS_proto.onwritestart =
	FS_proto.onprogress =
	FS_proto.onwrite =
	FS_proto.onabort =
	FS_proto.onerror =
	FS_proto.onwriteend =
	null;

	return saveAs;
}(
	typeof self !== "undefined" && self
	|| typeof window !== "undefined" && window
	|| this
));





	var my_icon =
	"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADMAAAAjCAMAAAAOqcKDAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAMAUExURQAAAAgUIQMiPgsiOggiPw0iPQ4kOh8hIhsiKhckPBUtNBMpPBkmOxwoPSUiHyAiIy0qLSMlPCEqOCsuMi0xNCkwPS0xOTErJjYuKTYuPTQ1OTc5PDg6Pjo8Pj4+PgAkQQAqRwAuSQAwSgczSxkwQic1SCw4Ris+TTM8Rzg6QD0+QTg/SRpBTCFITStOTj1HTz5ITkAoDkcyGEkyEUw7HUk4IUk4PlQ/P0A/Q1JGGFFIG1BIHVZIGlVJHVxGFFlEGV1HG1xIF1lIGE9JIEVGOE9OMU9PPFVIIlBNLlpILlZKP1hJPGFIGWBJI2FLKmFKLGRMKmpMImtMJmVNNG5QLWdQMGZUPWtQM2lSN3BSM3RXOnVYPkJBRENERUVBREZGR0JOT0dLT09GSE5IQE9PQk9PRUtLS0tLTElNT05NSE5OTkJNV0BRWEZSW0RSXlVMQFRPT15CQFNST1FVWVxTUVxYVlxdXUhddVBTbVZkbmBOTmVMSGZPTm1PRGNVTWdQTmpRTWZTUGdaUW9WUXNVTHhcRXlYSXBXUX1dUXtiTXtlU3hmWH9sW2RlY2Rmcm52e3BqZH5rZXRxZ3txaH53cWRwg2V2iYFfTIBeUYZhTYBtXopnU49sU4JwXoxzWpBpU5BrVJJsVZRvV5FwWJd1Wpt7XYJ1a455Y4p5b4d6cYl3cYt7cop+dox9dY1/eJl/ZYWAe4yBdo2Ce5aAbJCEfqOGZqeKZqqOaKuRbKySa6yUb62VcamYeK6bf7CbdLKhfoODg4KGi4SIjIqGg4iMkoaSm42QkJGHgpCIgpGKhpSKhpWMiJeRjJuRh5mRjZOVmJqWlIuXqYidrJOaopSfppujrJ2prZCptJqwt6iYjaiglaujmbSjgbKjibaqjbGkkbSql7askbevlratmLewl7ewmLixm6alpaimpqqoqKyqqq6urq60tqq0uLGqobSupLCvr7i0o7m3qbm3rbm4rrGxsbW0sbe3t7a5ubm3sLm4sbm4tbi4uP///wAAAAAAAHVRt0QAAAEAdFJOU////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////wBT9wclAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAAGXRFWHRTb2Z0d2FyZQBwYWludC5uZXQgNC4wLjE3M26fYwAAAuZJREFUSEvtU2VQlUEUvSZ2i92N3Tq2GIjdLbbyWaiI3a1gd6Nij90oPMVERewuEETEeAqi7+2d8dz9njP603+O4/nx7b1n79k9u3c/4j/HP6q5Oax8H89+Iw6aGXBuSIXBnv3cVyjEK3sbTd+aNI8dbNTXmsTh5ED7BD3z1ZEStY7lwF5EA/Zq3n6sMFF30agkjgIio3MiiPjijlTAEZPAt9Oa+MnliK5BY5+CmUo5PDyySslREJvKElXMP8pwSQkilLfgWydWNK8yI0yA5kwRosoX474+3miISHFULgx3lfq40yhENJtP1yDqvx8SdbIgURfcwefqqLikd96WHeFbDkpNVDdRKRW/PgVRW46YCroDCmyrEVyBJjg5UT05BXNw3zS5G4XzqfREblehUSfyOhVbxjZ/LFL7HXM0DNT8xGQPwP4LtIStS/0WPWAOKoXlqPmSN+pliOUWLJ3uQdQEBw0C3RH9sa1FsM/UOPB8QlIRCUKkP/xazHVl+wZID0NjnYX8sq79CVvAQC0QNAsRwj8dUcMEazWiWrBI1mmYwb6/4ssRr56pTA25YT3TXOhTnLwNps19dpm1zIt9fMz3c2O5t5EtVQnMtcDTiJKFuwXgI6egz2sQddd1bH9E5DWfE4Yakjw5NGNQMqLG4Q5zdUfCmrSWbP7QNDCfWfwcNO/wlz0oDJP8+wEPeBJzgTAnaCk0qeNys74SczQi1xjrPAzm7UdPlHMgiJqGF4CnultYws0WQLYDsfJF0Iptm9F9uiqzsoi2o81hQf1LENu2u5RB2tzbR1gKZ/txFwncx/iMLomxk5SxaU5bk38ucq6z5Cb8wEXOzOLIgCpoCBAtN+cpbrSGn03Pm9YsqDpe+q7OejsXNQnKeV+XsW0rzLnG6Fg0HLnOyJfRKWOG0rf1U2H1YqGRR4hMq95rAgj0MsY9NEOtYfUxzGK5cO+bTgTqw3WL5fydOEf6O0zNn+G/5m/WMP8AL8UQed9q5RYAAAAASUVORK5CYII=";


	mpld3.register_plugin("downloadprofile", DownloadProfile);
	DownloadProfile.prototype = Object.create(mpld3.Plugin.prototype);

	DownloadProfile.prototype.constructor = DownloadProfile;
	DownloadProfile.prototype.requiredProps = ["returnstr"];
	DownloadProfile.prototype.defaultProps = {};



	function DownloadProfile(fig, props){
		mpld3.Plugin.call(this, fig, props);
		var n = (this.props.returnstr).toString();



		var ResetButton = mpld3.ButtonFactory({
		buttonID: "reset",
		sticky: false,
		onActivate: function() {



			var blob = new Blob([n], {
				type: "text/plain;charset=utf-8;",
			});
			saveAs(blob, "Counts.csv");

		},
		icon: function(){return my_icon;
		}
		});
		this.fig.buttons.push(ResetButton);






	};

	"""

	def __init__(self,returnstr=2017,css=None):
		self.returnstr = returnstr
		self.css_ = css or ""

		self.dict_ = {"type": "downloadprofile",
					  "returnstr": returnstr}





class DownloadPNG(plugins.PluginBase):
	"""A Plugin to enable an HTML tooltip:
	formated text which hovers over points.

	Parameters
	----------
	points : matplotlib Collection or Line2D object
		The figure element to apply the tooltip to
	labels : list
		The labels for each point in points, as strings of unescaped HTML.
	hoffset, voffset : integer, optional
		The number of pixels to offset the tooltip text.  Default is
		hoffset = 0, voffset = 10
	css : str, optional
		css to be included, for styling the label html if desired
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> points = ax.plot(range(10), 'o')
	>>> labels = ['<h1>{title}</h1>'.format(title=i) for i in range(10)]
	>>> plugins.connect(fig, PointHTMLTooltip(points[0], labels))
	>>> fig_to_html(fig)
	"""

	JAVASCRIPT = """

	/* FileSaver.js
	 * A saveAs() FileSaver implementation.
	 * 1.3.8
	 * 2018-03-22 14:03:47
	 *
	 * By Eli Grey, https://eligrey.com
	 * License: MIT
	 *   See https://github.com/eligrey/FileSaver.js/blob/master/LICENSE.md
	 */

	/*global self */
	/*jslint bitwise: true, indent: 4, laxbreak: true, laxcomma: true, smarttabs: true, plusplus: true */

	/*! @source http://purl.eligrey.com/github/FileSaver.js/blob/master/src/FileSaver.js */

	var saveAs = saveAs || (function(view) {
	"use strict";
	// IE <10 is explicitly unsupported
	if (typeof view === "undefined" || typeof navigator !== "undefined" && /MSIE [1-9]\./.test(navigator.userAgent)) {
		return;
	}
	var
	doc = view.document
	// only get URL when necessary in case Blob.js hasn't overridden it yet
	, get_URL = function() {
		return view.URL || view.webkitURL || view;
	}
	, save_link = doc.createElementNS("http://www.w3.org/1999/xhtml", "a")
	, can_use_save_link = "download" in save_link
	, click = function(node) {
		var event = new MouseEvent("click");
		node.dispatchEvent(event);
	}
	, is_safari = /constructor/i.test(view.HTMLElement) || view.safari
	, is_chrome_ios =/CriOS\/[\d]+/.test(navigator.userAgent)
	, setImmediate = view.setImmediate || view.setTimeout
	, throw_outside = function(ex) {
		setImmediate(function() {
			throw ex;
		}, 0);
	}
	, force_saveable_type = "application/octet-stream"
	// the Blob API is fundamentally broken as there is no "downloadfinished" event to subscribe to
	, arbitrary_revoke_timeout = 1000 * 40 // in ms
	, revoke = function(file) {
		var revoker = function() {
			if (typeof file === "string") { // file is an object URL
				get_URL().revokeObjectURL(file);
			} else { // file is a File
				file.remove();
			}
		};
		setTimeout(revoker, arbitrary_revoke_timeout);
	}
	, dispatch = function(filesaver, event_types, event) {
		event_types = [].concat(event_types);
		var i = event_types.length;
		while (i--) {
			var listener = filesaver["on" + event_types[i]];
			if (typeof listener === "function") {
				try {
					listener.call(filesaver, event || filesaver);
				} catch (ex) {
					throw_outside(ex);
				}
			}
		}
	}
	, auto_bom = function(blob) {
		// prepend BOM for UTF-8 XML and text/* types (including HTML)
		// note: your browser will automatically convert UTF-16 U+FEFF to EF BB BF
		if (/^\s*(?:text\/\S*|application\/xml|\S*\/\S*\+xml)\s*;.*charset\s*=\s*utf-8/i.test(blob.type)) {
			//return new Blob([String.fromCharCode(0xFEFF), blob], {type: blob.type});
			return new Blob(["\uFEFF", blob], {type: blob.type});
		}
		return blob;
	}
	, FileSaver = function(blob, name, no_auto_bom) {
		if (!no_auto_bom) {
			blob = auto_bom(blob);
		}
		// First try a.download, then web filesystem, then object URLs
		var
		filesaver = this
		, type = blob.type
		, force = type === force_saveable_type
		, object_url
		, dispatch_all = function() {
			dispatch(filesaver, "writestart progress write writeend".split(" "));
		}
		// on any filesys errors revert to saving with object URLs
		, fs_error = function() {
			if ((is_chrome_ios || (force && is_safari)) && view.FileReader) {
				// Safari doesn't allow downloading of blob urls
				var reader = new FileReader();
				reader.onloadend = function() {
					var url = is_chrome_ios ? reader.result : reader.result.replace(/^data:[^;]*;/, 'data:attachment/file;');
					var popup = view.open(url, '_blank');
					if(!popup) view.location.href = url;
							   url=undefined; // release reference before dispatching
							   filesaver.readyState = filesaver.DONE;
					dispatch_all();
				};
				reader.readAsDataURL(blob);
				filesaver.readyState = filesaver.INIT;
				return;
			}
			// don't create more object URLs than needed
			if (!object_url) {
				object_url = get_URL().createObjectURL(blob);
			}
			if (force) {
				view.location.href = object_url;
			} else {
				var opened = view.open(object_url, "_blank");
				if (!opened) {
					// Apple does not allow window.open, see https://developer.apple.com/library/safari/documentation/Tools/Conceptual/SafariExtensionGuide/WorkingwithWindowsandTabs/WorkingwithWindowsandTabs.html
					view.location.href = object_url;
				}
			}
			filesaver.readyState = filesaver.DONE;
			dispatch_all();
			revoke(object_url);
		}
		;
		filesaver.readyState = filesaver.INIT;

		if (can_use_save_link) {
			object_url = get_URL().createObjectURL(blob);
			setImmediate(function() {
				save_link.href = object_url;
				save_link.download = name;
				click(save_link);
				dispatch_all();
				revoke(object_url);
				filesaver.readyState = filesaver.DONE;
			}, 0);
			return;
		}

		fs_error();
	}
	, FS_proto = FileSaver.prototype
	, saveAs = function(blob, name, no_auto_bom) {
		return new FileSaver(blob, name || blob.name || "download", no_auto_bom);
	}
	;

	// IE 10+ (native saveAs)
	if (typeof navigator !== "undefined" && navigator.msSaveOrOpenBlob) {
		return function(blob, name, no_auto_bom) {
			name = name || blob.name || "download";

			if (!no_auto_bom) {
				blob = auto_bom(blob);
			}
			return navigator.msSaveOrOpenBlob(blob, name);
		};
	}

	// todo: detect chrome extensions & packaged apps
	//save_link.target = "_blank";

	FS_proto.abort = function(){};
	FS_proto.readyState = FS_proto.INIT = 0;
	FS_proto.WRITING = 1;
	FS_proto.DONE = 2;

	FS_proto.error =
	FS_proto.onwritestart =
	FS_proto.onprogress =
	FS_proto.onwrite =
	FS_proto.onabort =
	FS_proto.onerror =
	FS_proto.onwriteend =
	null;

	return saveAs;
	}(
		typeof self !== "undefined" && self
		|| typeof window !== "undefined" && window
		|| this
	));



	function getSVGString( svgNode ) {
		svgNode.setAttribute('xlink', 'http://www.w3.org/1999/xlink');
		var cssStyleText = getCSSStyles( svgNode );
		appendCSS( cssStyleText, svgNode );
	
		var serializer = new XMLSerializer();
		var svgString = serializer.serializeToString(svgNode);
		svgString = svgString.replace(/(\w+)?:?xlink=/g, 'xmlns:xlink='); // Fix root xlink without namespace
		svgString = svgString.replace(/NS\d+:href/g, 'xlink:href'); // Safari NS namespace fix
	
		return svgString;
	
		function getCSSStyles( parentElement ) {
			var selectorTextArr = [];
	
			// Add Parent element Id and Classes to the list
			selectorTextArr.push( '#'+parentElement.id );
			for (var c = 0; c < parentElement.classList.length; c++)
					if ( !contains('.'+parentElement.classList[c], selectorTextArr) )
						selectorTextArr.push( '.'+parentElement.classList[c] );
	
			// Add Children element Ids and Classes to the list
			var nodes = parentElement.getElementsByTagName("*");
			for (var i = 0; i < nodes.length; i++) {
				var id = nodes[i].id;
				if ( !contains('#'+id, selectorTextArr) )
					selectorTextArr.push( '#'+id );
	
				var classes = nodes[i].classList;
				for (var c = 0; c < classes.length; c++)
					if ( !contains('.'+classes[c], selectorTextArr) )
						selectorTextArr.push( '.'+classes[c] );
			}
	
			// Extract CSS Rules
			var extractedCSSText = "";
			for (var i = 0; i < document.styleSheets.length; i++) {
				var s = document.styleSheets[i];
				
				try {
					if(!s.cssRules) continue;
				} catch( e ) {
						if(e.name !== 'SecurityError') throw e; // for Firefox
						continue;
					}
	
				var cssRules = s.cssRules;
				for (var r = 0; r < cssRules.length; r++) {
					if ( contains( cssRules[r].selectorText, selectorTextArr ) )
						extractedCSSText += cssRules[r].cssText;
				}
			}
			
	
			return extractedCSSText;
	
			function contains(str,arr) {
				return arr.indexOf( str ) === -1 ? false : true;
			}
	
		}
	
		function appendCSS( cssText, element ) {
			var styleElement = document.createElement("style");
			styleElement.setAttribute("type","text/css"); 
			styleElement.innerHTML = cssText;
			var refNode = element.hasChildNodes() ? element.children[0] : null;
			element.insertBefore( styleElement, refNode );
		}
	}

	function svgString2Image( svgString, width, height, format, callback ) {
		var format = format ? format : 'png';
	
		var imgsrc = 'data:image/svg+xml;base64,'+ btoa( unescape( encodeURIComponent( svgString ) ) ); // Convert SVG string to data URL
	
		var canvas = document.createElement("canvas");
		var context = canvas.getContext("2d");
	
		canvas.width = width;
		canvas.height = height;
	
		var image = new Image();
		image.onload = function() {
			context.clearRect ( 0, 0, width, height );
			context.drawImage(image, 0, 0, width, height);
	
			canvas.toBlob( function(blob) {
				var filesize = Math.round( blob.length/1024 ) + ' KB';
				if ( callback ) callback( blob, filesize );
			});
	
			
		};
	
		image.src = imgsrc;
	}


	var my_save_icon = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADMAAAAjCAMAAAAOqcKDAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAMAUExURQAAAP///wAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGd27GMAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAZdEVYdFNvZnR3YXJlAHBhaW50Lm5ldCA0LjAuMTczbp9jAAAAk0lEQVRIS+2LSQ6AMAwD6f8/DU7iLFWpCkiIAz40Ezeztev5nTedjbGiotIR6XVPjhWE4kir6RxcOGRHOsYc78dOVMDkGOAVCMcL5KBVB5MZOni8iAPmmWNTB59odCLmMDcc/uhIeOqgnzmIsP8wXDETYhJXHTsROnVMCnSaON1lZM2p1sypP1gQoDaX8jtfdlrbAc3PBRcmkelCAAAAAElFTkSuQmCC";

	mpld3.register_plugin("downloadpng", DownloadPNG);
	DownloadPNG.prototype = Object.create(mpld3.Plugin.prototype);

	DownloadPNG.prototype.constructor = DownloadPNG;
	DownloadPNG.prototype.requiredProps = ["returnstr"];
	DownloadPNG.prototype.defaultProps = {};


	function saveCanvas(x_canvas){
		x_canvas.toBlob(function(blob) {
			saveAs(
				blob
				, "screenshot.png"
			);
		}, "image/png");
	}
	function DownloadPNG(fig, props){
		mpld3.Plugin.call(this, fig, props);
		var n = (this.props.returnstr).toString();

		var ResetButton = mpld3.ButtonFactory({
		buttonID: "reset",
		sticky: false,
		onActivate: function() {
			//var blob = new Blob([n], {
			//	type: "text/plain;charset=utf-8;",
			//});
			//saveAs(blob, "Image.html");
			//base_image = new Image();
			//base_image.src = 'data:image/png;base64,'+n
			//var canvas = document.createElement('canvas');
			//canvas.id     = "YourCanvas";
			//var canvas = document.getElementById('YourCanvas');
			//context = canvas.getContext('2d');
			// Draw image within
			//context.drawImage(base_image, 0,0);
			//canvas.toBlob(function(blob) {
			//	saveAs(blob, "Image.png");
			//}, "image/png");
			
			var svgs = document.getElementsByClassName("mpld3-figure");
			svgString = getSVGString(svgs[0])
			svgString2Image( svgString, 2300, 1200, 'png', save ); // passes Blob and filesize String to the callback
			
			function save( dataBlob, filesize ){
				saveAs( dataBlob, n ); // FileSaver.js function
			}
			console.log("svg string"+svgString);
			
			
			base_image = new Image();
			base_image.src ='data:image/png;base64,'+n;
			base_image.onload = function(){
				var canvas = document.createElement('canvas');
				canvas.width = 2300;
				canvas.height = 1100;
				context = canvas.getContext('2d');
				// Draw image within
				context.drawImage(base_image, 0,0);
				// Save the canvas
				//saveCanvas(canvas);
				// Remove it
				
			};
		},
		icon: function(){return my_save_icon;}
		});
		this.fig.buttons.push(ResetButton);
	};

	"""

	def __init__(self,returnstr=2017,css=None):
		self.returnstr = returnstr
		self.css_ = css or ""

		self.dict_ = {"type": "downloadpng",
					  "returnstr": returnstr}
