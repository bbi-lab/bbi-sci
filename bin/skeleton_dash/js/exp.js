var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

function Header(props) {
  return React.createElement(
    "nav",
    { className: "navbar navbar-expand-md sticky-top navbar-light", style: { backgroundColor: "#e3f2fd" } },
    React.createElement(
      "div",
      { className: "navbar-collapse collapse w-100 order-1 order-md-0 dual-collapse2" },
      React.createElement(
        "ul",
        { className: "navbar-nav mr-auto" },
        React.createElement("img", { src: "img/bbi_icon.png", height: "70", className: "d-inline-block align-top", alt: "" })
      )
    ),
    React.createElement(
      "div",
      { className: "mx-auto order-0" },
      React.createElement(
        "a",
        { className: "navbar-brand mx-auto", href: "#" },
        "Experiment ",
        props.run_name,
        " QC Dashboard"
      )
    ),
    React.createElement("div", { className: "navbar-collapse collapse w-100 order-3 dual-collapse2" })
  );
}

function Sample(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: props.sample_id, role: "tabpanel",
      "aria-labelledby": props.sample_id },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "lig-name" },
        props.sample_id
      )
    ),
    React.createElement(
      "nav",
      null,
      React.createElement(
        "div",
        { className: "nav nav-tabs", id: "nav" + props.sample_id + "-tab", role: "tablist" },
        props.sample_id == "Barnyard" && React.createElement(
          "a",
          {
            className: "nav-item nav-link active",
            id: "nav" + props.sample_id + "-barn-tab",
            "data-toggle": "tab", href: "#nav" + props.sample_id + "-barn",
            role: "tab", "aria-controls": "nav" + props.sample_id + "-barn",
            "aria-selected": "true" },
          "Barnyard"
        ),
        props.sample_id == "Barnyard" ? React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + props.sample_id + "-knee-tab",
            "data-toggle": "tab", href: "#nav" + props.sample_id + "-knee",
            role: "tab", "aria-controls": "nav" + props.sample_id + "-knee",
            "aria-selected": "true" },
          "Knee plot"
        ) : React.createElement(
          "a",
          {
            className: "nav-item nav-link active",
            id: "nav" + props.sample_id + "-knee-tab",
            "data-toggle": "tab", href: "#nav" + props.sample_id + "-knee",
            role: "tab", "aria-controls": "nav" + props.sample_id + "-knee",
            "aria-selected": "true" },
          "Knee plot"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + props.sample_id + "-scrub-tab",
            "data-toggle": "tab", href: "#nav" + props.sample_id + "-scrub",
            role: "tab", "aria-controls": "nav" + props.sample_id + "-scrub",
            "aria-selected": "false" },
          "Scrublet"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + props.sample_id + "-cellqc-tab",
            "data-toggle": "tab", href: "#nav" + props.sample_id + "-cellqc",
            role: "tab", "aria-controls": "nav" + props.sample_id + "-cellqc",
            "aria-selected": "false" },
          "Cell QC"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + props.sample_id + "-umap-tab",
            "data-toggle": "tab", href: "#nav" + props.sample_id + "-umap",
            role: "tab", "aria-controls": "nav" + props.sample_id + "-umap",
            "aria-selected": "false" },
          "UMAP"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + props.sample_id + "-stats-tab",
            "data-toggle": "tab", href: "#nav" + props.sample_id + "-stats",
            role: "tab", "aria-controls": "nav" + props.sample_id + "-stats",
            "aria-selected": "false" },
          "Sample Stats"
        )
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-tabContent" },
      props.sample_id == "Barnyard" && React.createElement(BarnyardPane, { sample_id: props.sample_id, className: "tab-pane fade show active" }),
      props.sample_id == "Barnyard" ? React.createElement(KneePane, { sample_id: props.sample_id, className: "tab-pane fade" }) : React.createElement(KneePane, { sample_id: props.sample_id, className: "tab-pane fade show active" }),
      React.createElement(ScrubPane, { sample_id: props.sample_id, sample_stats: run_data.sample_stats }),
      React.createElement(QCPane, { sample_id: props.sample_id }),
      React.createElement(UMAPPane, { sample_id: props.sample_id }),
      React.createElement(StatsPane, { sample_id: props.sample_id, sample_stats: run_data.sample_stats })
    )
  );
}

function Pane(props) {
  return React.createElement(
    "div",
    { className: props.className, id: props.id, role: "tabpanel", "aria-labelledby": props.tag },
    React.createElement(
      "p",
      null,
      props.text.map(function (text, index) {
        return typeof text == "string" ? React.createElement(
          "span",
          { key: index },
          text
        ) : React.createElement(
          "a",
          { key: index, href: text.link },
          text.label
        );
      })
    ),
    React.createElement("img", { src: props.plot, className: "rounded mx-auto d-block", alt: "...", style: { maxHeight: "50vh", width: "auto" } })
  );
}

function TitleRow(props) {
  return React.createElement(
    "th",
    { scope: "col" },
    props.samp
  );
}

function RegRow(props) {
  return React.createElement(
    "td",
    null,
    props.val
  );
}

function StatsPane(props) {
  var sample_stat = props.sample_stats[props.sample_id];
  var stats_list = ["Total Reads", "Total UMIs", "Duplication Rate", "Cells with >100 UMIs", "Cells with >1000 UMIs"];
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "nav" + props.sample_id + "-stats", role: "tabpanel", "aria-labelledby": "nav" + props.sample_id + "-stats-tab" },
    React.createElement(
      "table",
      { className: "table table-hover" },
      React.createElement(
        "thead",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement("th", { scope: "col" }),
          stats_list.map(function (item, index) {
            return React.createElement(TitleRow, { key: index, samp: item });
          })
        )
      ),
      React.createElement(
        "tbody",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            props.sample_id
          ),
          React.createElement(RegRow, { val: sample_stat.Total_reads }),
          React.createElement(RegRow, { val: sample_stat.Total_UMIs }),
          React.createElement(RegRow, { val: sample_stat.Duplication_rate }),
          React.createElement(RegRow, { val: sample_stat.Cells_100_UMIs }),
          React.createElement(RegRow, { val: sample_stat.Cells_1000_UMIs })
        )
      )
    )
  );
}

function BarnyardPane(props) {
  return React.createElement(Pane, {
    className: props.className,
    id: "nav" + props.sample_id + "-barn",
    tag: "nav" + props.sample_id + "-barn-tab",
    text: ["Collision rate: " + run_data.barn_collision],
    plot: "img/Barnyard_plot.png"
  });
}
function QCPane(props) {
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + props.sample_id + "-cellqc",
    tag: "nav" + props.sample_id + "-cellqc-tab",
    text: [''],
    plot: "img/" + props.sample_id + "_cell_qc.png"
  });
}
function ScrubPane(props) {
  var sample_stat = props.sample_stats[props.sample_id];
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + props.sample_id + "-scrub",
    tag: "nav" + props.sample_id + "-scrub-tab",
    text: ['Doublet count: ' + sample_stat.Doublet_Number + "\n\nDoublet rate: " + sample_stat.Doublet_Percent],
    plot: "img/" + props.sample_id + "_scrublet_hist.png"
  });
}
function KneePane(props) {
  return React.createElement(Pane, {
    className: props.className,
    id: "nav" + props.sample_id + "-knee",
    tag: "nav" + props.sample_id + "-knee-tab",
    text: [''],
    plot: "img/" + props.sample_id + "_knee_plot.png"
  });
}
function UMAPPane(props) {
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + props.sample_id + "-umap",
    tag: "nav" + props.sample_id + "-umap-tab",
    text: [''],
    plot: "img/" + props.sample_id + "_UMAP.png"
  });
}

function SamplePill(props) {
  return React.createElement(
    "a",
    { className: "nav-link", id: props.sample_id + "-tab", "data-toggle": "pill", href: "#" + props.sample_id, role: "tab",
      "aria-controls": props.sample_id, "aria-selected": "false" },
    props.sample_id
  );
}

//https://www.florin-pop.com/blog/2019/07/sort-table-data-with-react/

var tableData = Object.values(run_data['sample_stats']);
console.log(tableData);

var sortTypes = {
  reads_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return a.Total_reads - b.Total_reads;
    }
  },
  reads_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return b.Total_reads - a.Total_reads;
    }
  },
  umis_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return a.Total_UMIs - b.Total_UMIs;
    }
  },
  umis_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return b.Total_UMIs - a.Total_UMIs;
    }
  },
  sample_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return ('' + a.Sample).localeCompare(b.Sample);
    }
  },
  sample_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return ('' + b.Sample).localeCompare(a.Sample);
    }
  },
  dup_rate_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return parseFloat(a.Duplication_rate) - parseFloat(b.Duplication_rate);
    }
  },
  dup_rate_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return parseFloat(b.Duplication_rate) - parseFloat(a.Duplication_rate);
    }
  },
  doub_rate_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return a.Doublet_Percent == "Fail" ? -1 : parseFloat(a.Doublet_Percent) - parseFloat(b.Doublet_Percent);
    }
  },
  doub_rate_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return b.Doublet_Percent == "Fail" ? -1 : parseFloat(b.Doublet_Percent) - parseFloat(a.Doublet_Percent);
    }
  },
  c100_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return b.Cells_100_UMIs - a.Cells_100_UMIs;
    }
  },
  c100_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return a.Cells_100_UMIs - b.Cells_100_UMIs;
    }
  },
  c1000_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return b.Cells_1000_UMIs - a.Cells_1000_UMIs;
    }
  },
  c1000_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return a.Cells_1000_UMIs - b.Cells_1000_UMIs;
    }
  },
  default: {
    class: 'sort',
    fn: function fn(a, b) {
      return a;
    }
  }
};

var Table = function (_React$Component) {
  _inherits(Table, _React$Component);

  function Table() {
    var _ref;

    var _temp, _this, _ret;

    _classCallCheck(this, Table);

    for (var _len = arguments.length, args = Array(_len), _key = 0; _key < _len; _key++) {
      args[_key] = arguments[_key];
    }

    return _ret = (_temp = (_this = _possibleConstructorReturn(this, (_ref = Table.__proto__ || Object.getPrototypeOf(Table)).call.apply(_ref, [this].concat(args))), _this), _this.state = {
      currentSort: 'sample_up'
    }, _this.onSortReads = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'reads_down') nextSort = 'reads_up';else if (currentSort === 'reads_up') nextSort = 'reads_down';else nextSort = 'reads_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortSample = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'sample_down') nextSort = 'sample_up';else if (currentSort === 'sample_up') nextSort = 'sample_down';else nextSort = 'sample_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortDoub = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'doub_rate_down') nextSort = 'doub_rate_up';else if (currentSort === 'doub_rate_up') nextSort = 'doub_rate_down';else nextSort = 'doub_rate_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortUMIs = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'umis_down') nextSort = 'umis_up';else if (currentSort === 'umis_up') nextSort = 'umis_down';else nextSort = 'umis_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortDup = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'dup_rate_down') nextSort = 'dup_rate_up';else if (currentSort === 'dup_rate_up') nextSort = 'dup_rate_down';else nextSort = 'dup_rate_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortc100 = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'c100_down') nextSort = 'c100_up';else if (currentSort === 'c100_up') nextSort = 'c100_down';else nextSort = 'c100_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortc1000 = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'c1000_down') nextSort = 'c1000_up';else if (currentSort === 'c1000_up') nextSort = 'c1000_down';else nextSort = 'c1000_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _temp), _possibleConstructorReturn(_this, _ret);
  }

  // declaring the default state


  // method called every time the sort button is clicked
  // it will change the currentSort value to the next one


  _createClass(Table, [{
    key: "render",
    value: function render() {
      var data = this.props.data;
      var currentSort = this.state.currentSort;

      return React.createElement(
        "div",
        { className: "tab-pane fade show active", id: "summary", role: "tabpanel" },
        React.createElement(
          "div",
          { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
          React.createElement(
            "h1",
            { className: "h3", id: "lig-name" },
            "Summary Table"
          )
        ),
        data.length > 0 && React.createElement(
          "table",
          { className: "table table-hover table-responsive summary-table" },
          React.createElement(
            "thead",
            null,
            React.createElement(
              "tr",
              null,
              React.createElement(
                "th",
                null,
                "Sample",
                React.createElement(
                  "button",
                  { onClick: this.onSortSample, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Total reads",
                React.createElement(
                  "button",
                  { onClick: this.onSortReads, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Total UMIs",
                React.createElement(
                  "button",
                  { onClick: this.onSortUMIs, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Duplication rate",
                React.createElement(
                  "button",
                  { onClick: this.onSortDup, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Doublet rate",
                React.createElement(
                  "button",
                  { onClick: this.onSortDoub, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Cells with >100 UMIs",
                React.createElement(
                  "button",
                  { onClick: this.onSortc100, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Cells with >1000 UMIs",
                React.createElement(
                  "button",
                  { onClick: this.onSortc1000, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              )
            )
          ),
          React.createElement(
            "tbody",
            null,
            [].concat(_toConsumableArray(data)).sort(sortTypes[currentSort].fn).map(function (p) {
              return React.createElement(
                "tr",
                null,
                React.createElement(
                  "td",
                  null,
                  p.Sample
                ),
                React.createElement(
                  "td",
                  null,
                  p.Total_reads
                ),
                React.createElement(
                  "td",
                  null,
                  p.Total_UMIs
                ),
                React.createElement(
                  "td",
                  null,
                  p.Duplication_rate
                ),
                React.createElement(
                  "td",
                  null,
                  p.Doublet_Percent
                ),
                React.createElement(
                  "td",
                  null,
                  p.Cells_100_UMIs
                ),
                React.createElement(
                  "td",
                  null,
                  p.Cells_1000_UMIs
                )
              );
            })
          )
        )
      );
    }
  }]);

  return Table;
}(React.Component);

function ExperimentPage(props) {
  return React.createElement(
    "span",
    null,
    React.createElement(Header, { run_name: props.run_name }),
    React.createElement(
      "div",
      { className: "container-fluid" },
      React.createElement(
        "div",
        { className: "row" },
        React.createElement(
          "nav",
          { className: "col-md-2 d-none d-md-block bg-light sidebar" },
          React.createElement(
            "div",
            { className: "sidebar-sticky" },
            React.createElement(
              "div",
              { className: "nav flex-column nav-pills", id: "v-pills-tab", role: "tablist", "aria-orientation": "vertical" },
              React.createElement(
                "a",
                { className: "nav-link active", id: "summary-tab", "data-toggle": "pill", href: "#summary", role: "tab", "aria-controls": "summary", "aria-selected": "true" },
                "Summary Table"
              ),
              props.samp_pills
            )
          )
        ),
        React.createElement(
          "main",
          { role: "main", className: "col-md-9 ml-sm-auto col-lg-10 px-4", style: { paddingTop: "15px" } },
          React.createElement(
            "div",
            { className: "tab-content", id: "nav-tabContent" },
            React.createElement(Table, { data: tableData }),
            props.samp_list
          )
        )
      )
    )
  );
}

var sampList = run_data.sample_list.map(function (samp) {
  return React.createElement(Sample, { key: samp, sample_id: samp });
});

var sampPills = run_data.sample_list.map(function (samp) {
  return React.createElement(SamplePill, { key: samp, sample_id: samp });
});

ReactDOM.render(React.createElement(ExperimentPage, { samp_list: sampList, samp_pills: sampPills, run_name: run_data.run_name }), document.getElementById('exp_page'));