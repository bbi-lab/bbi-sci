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
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: safe_name, role: "tabpanel",
      "aria-labelledby": safe_name },
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
        { className: "nav nav-tabs", id: "nav" + safe_name + "-tab", role: "tablist" },
        props.sample_id == "Barnyard" && React.createElement(
          "a",
          {
            className: "nav-item nav-link active",
            id: "nav" + safe_name + "-barn-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-barn",
            role: "tab", "aria-controls": "nav" + safe_name + "-barn",
            "aria-selected": "true" },
          "Barnyard"
        ),
        props.sample_id == "Barnyard" ? React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-knee-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-knee",
            role: "tab", "aria-controls": "nav" + safe_name + "-knee",
            "aria-selected": "true" },
          "Knee plot"
        ) : React.createElement(
          "a",
          {
            className: "nav-item nav-link active",
            id: "nav" + safe_name + "-knee-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-knee",
            role: "tab", "aria-controls": "nav" + safe_name + "-knee",
            "aria-selected": "true" },
          "Knee plot"
        ),

        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-umis-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-umis",
            role: "tab", "aria-controls": "nav" + safe_name + "-umis",
            "aria-selected": "false" },
          "UMI Plots"
        ),

        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-genes-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-genes",
            role: "tab", "aria-controls": "nav" + safe_name + "-genes",
            "aria-selected": "false" },
          "Genes By UMIs"
        ),

        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-pseudobulk-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-pseudobulk",
            role: "tab", "aria-controls": "nav" + safe_name + "-pseudobulk",
            "aria-selected": "false" },
          "Barcode Correlations Heatmap"
        ),

        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-pseudobulk-hist-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-pseudobulk-hist",
            role: "tab", "aria-controls": "nav" + safe_name + "-pseudobulk-hist",
            "aria-selected": "false" },
          "Barcode Correlations Histogram"
        ),


        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-hash-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-hash",
            role: "tab", "aria-controls": "nav" + safe_name + "-hash",
            "aria-selected": "false" },
          "Hash Plots"
        ),
        props.garnett_model != null && React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-garnett-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-garnett",
            role: "tab", "aria-controls": "nav" + safe_name + "-garnett",
            "aria-selected": "false" },
          "Garnett"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-stats-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-stats",
            role: "tab", "aria-controls": "nav" + safe_name + "-stats",
            "aria-selected": "false" },
          "Sample Stats"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-readmetrics-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-readmetrics",
            role: "tab", "aria-controls": "nav" + safe_name + "-readmetrics",
            "aria-selected": "false" },
          "Read Metrics"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-fulllog-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-fulllog",
            role: "tab", "aria-controls": "nav" + safe_name + "-fulllog",
            "aria-selected": "false" },
          "Full Log"
        )
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-tabContent" },
      props.sample_id == "Barnyard" && React.createElement(BarnyardPane, { sample_id: props.sample_id, className: "tab-pane fade show active" }),
      props.sample_id == "Barnyard" ? React.createElement(KneePane, { sample_id: props.sample_id, className: "tab-pane fade" }) : React.createElement(KneePane, { sample_id: props.sample_id, className: "tab-pane fade show active" }),
      props.garnett_model != null && React.createElement(GarnettPane, { sample_id: props.sample_id, garnett_model: props.garnett_model }),
      React.createElement(StatsPane, { sample_id: props.sample_id, sample_stats: run_data.sample_stats }),
      React.createElement(UMIPane, { sample_id: props.sample_id }),
      React.createElement(GeneByUMIPane, { sample_id: props.sample_id }),
      React.createElement(PseudobulkPane, { sample_id: props.sample_id }),
      React.createElement(PseudobulkHistPane, { sample_id: props.sample_id }),
      React.createElement(HashPane, { sample_id: props.sample_id }),
      React.createElement(ReadMetricsPane, { sample_id: props.sample_id, log: log_data[props.sample_id] }),
      React.createElement(FullLogPane, { sample_id: props.sample_id, log: full_log_data[props.sample_id] }),
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
  // var stats_list = ["Total Reads", "Total UMIs", "Median UMIs", "Median Mitochondrial UMIs", "Duplication Rate", ",Cells with >100 UMIs", "Cells with >1000 UMIs", "Cells with FDR<=.01", "Cells with FDR<=.001"];
  var stats_list = ["Total Reads", "Total UMIs", "Median UMIs", "Median Mitochondrial UMIs", "Duplication Rate", ",Cells with >100 UMIs", "Cells with FDR<=.01"];


  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "nav" + safe_name + "-stats", role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-stats-tab" },
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
          React.createElement(RegRow, { val: sample_stat.Median_UMIs }),
          React.createElement(RegRow, { val: sample_stat.Median_Mitochondrial_UMIs_Percent}),
          React.createElement(RegRow, { val: sample_stat.Duplication_rate }),
          React.createElement(RegRow, { val: sample_stat.Cells_100_UMIs }),
          // React.createElement(RegRow, { val: sample_stat.Cells_1000_UMIs }),
          React.createElement(RegRow, { val: sample_stat.Cells_FDR_p01 })
          // React.createElement(RegRow, { val: sample_stat.Cells_FDR_p001 })
        )
      )
    )
  );
}

function CodeChunk(props) {
  return React.createElement(
    "pre",
    { style: { paddingLeft: '20px', paddingRight: '20px' } },
    React.createElement(
      "code",
      null,
      '\n' + props.text + '\n\n'
    )
  );
}

function ReadMetricsPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "nav" + safe_name + "-readmetrics", role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-readmetrics-tab" },
    React.createElement(CodeChunk, { text: props.log })
  );
}

function FullLogPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "nav" + safe_name + "-fulllog", role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-fulllog-tab" },
    React.createElement(CodeChunk, { text: props.log })
  );
}

function BarnyardPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement(Pane, {
    className: props.className,
    id: "nav" + safe_name + "-barn",
    tag: "nav" + safe_name + "-barn-tab",
    text: ["Collision rate: " + run_data.barn_collision],
    plot: "img/Barnyard_plot.png"
  });
}
// function QCPane(props) {
//   var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
//   return React.createElement(Pane, {
//     className: "tab-pane fade",
//     id: "nav" + safe_name + "-cellqc",
//     tag: "nav" + safe_name + "-cellqc-tab",
//     text: [''],
//     plot: "img/" + props.sample_id + "_cell_qc.png"
//   });
// }
// function ScrubPane(props) {
//   var sample_stat = props.sample_stats[props.sample_id];
//   var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
//   return React.createElement(Pane, {
//     className: "tab-pane fade",
//     id: "nav" + safe_name + "-scrub",
//     tag: "nav" + safe_name + "-scrub-tab",
//     // text: ['Doublet count: ' + sample_stat.Doublet_Number + "\n\nDoublet rate: " + sample_stat.Doublet_Percent],
//     text: [],
//     plot: "img/" + props.sample_id + "_scrublet_hist.png"
//   });
// }
function KneePane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement(Pane, {
    className: props.className,
    id: "nav" + safe_name + "-knee",
    tag: "nav" + safe_name + "-knee-tab",
    text: [''],
    plot: "img/" + props.sample_id + "_knee_plot.png"
  });
}

// function UMAPPane(props) {
//   var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
//   return React.createElement(Pane, {
//     className: "tab-pane fade",
//     id: "nav" + safe_name + "-umap",
//     tag: "nav" + safe_name + "-umap-tab",
//     text: [''],
//     plot: "img/" + props.sample_id + "_UMAP.png"
//   });
// }


function GarnettPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "nav" + safe_name + "-garnett", role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-garnett-tab" },
    Array.isArray(props.garnett_model) ? props.garnett_model.map(function (model, index) {
      return React.createElement(
        "span",
        { key: index },
        React.createElement(
          "p",
          null,
          'Garnett model run: ' + model
        ),
        React.createElement("img", { src: "img/" + props.sample_id + "_" + model + "_Garnett.png", className: "rounded mx-auto d-block", alt: "...", style: { maxHeight: "50vh", width: "auto" } })
      );
    }) : props.garnett_model == "no_cells" ? React.createElement(
      "p",
      null,
      "No cells to apply Garnett models."
    ) : React.createElement(
      "span",
      null,
      React.createElement(
        "p",
        null,
        'Garnett model run: ' + props.garnett_model
      ),
      React.createElement("img", { src: "img/" + props.sample_id + "_" + props.garnett_model + "_Garnett.png", className: "rounded mx-auto d-block", alt: "...", style: { maxHeight: "50vh", width: "auto" } })
    )
  );
}
function UMIPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement("div", { 
    className: "tab-pane fade",
    id: "nav" + safe_name + "-umis",
    role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-umis-tab" },

    React.createElement("img", { 
      src: "img/" + props.sample_id + "_umi.png", 
      className: "rounded mx-auto d-block", 
      alt: "...",
      style: { maxHeight: "100vh", width: "auto" } 
    }),
  );
}
function GeneByUMIPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement("div", { 
    className: "tab-pane fade",
    id: "nav" + safe_name + "-genes",
    role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-genes-tab" },

    React.createElement("img", { 
      src: "img/" + props.sample_id + "_genes_by_umi.png", 
      className: "rounded mx-auto d-block", 
      alt: "...",
      style: { maxHeight: "100vh", width: "auto" } 
    }),
  );
}

function PseudobulkPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");

  return React.createElement("div", { 
    className: "tab-pane fade",
    id: "nav" + safe_name + "-pseudobulk",
    role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-pseudobulk-tab" },

    React.createElement("img", { 
      src: "img/" + props.sample_id + "_pseudobulk_heatmap.png", 
      className: "rounded mx-auto d-block", 
      alt: "...",
      style: { maxHeight: "200vh", width: "auto" } 
    }),
  );
}

function PseudobulkHistPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");

  return React.createElement("div", { 
    className: "tab-pane fade",
    id: "nav" + safe_name + "-pseudobulk-hist",
    role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-pseudobulk-hist-tab" },

    React.createElement("img", { 
      src: "img/" + props.sample_id + "_pseudobulk_histogram.png", 
      className: "rounded mx-auto d-block", 
      alt: "...",
      style: { maxHeight: "120vh", width: "auto" } 
    }),
  );
}

function HashPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement("div", { 
    className: "tab-pane fade",
    id: "nav" + safe_name + "-hash",
    role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-hash-tab" },

    React.createElement("img", { 
      src: "img/" + props.sample_id + "_hash_plots.png", 
      className: "rounded mx-auto d-block", 
      alt: "...",
      style: { maxHeight: "100vh", width: "auto" } 
    }),
  );
}


// function WellCheckPane(props) {
//   var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");

//   return React.createElement("div", { 
//     className: "tab-pane fade",
//     id: "nav" + safe_name + "-wellcheck",
//     role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-wellcheck-tab" },

//     React.createElement("img", { 
//       src: "img/" + props.sample_id + "_wellcheck.png", 
//       className: "rounded mx-auto d-block", 
//       alt: "...",
//       style: { maxHeight: "100vh", width: "auto" } 
//     }),
//   );
// }




function SamplePill(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "");
  return React.createElement(
    "a",
    { className: "nav-link", id: safe_name + "-tab", "data-toggle": "pill", href: "#" + safe_name, role: "tab",
      "aria-controls": safe_name, "aria-selected": "false" },
    props.sample_id
  );
}

//https://www.florin-pop.com/blog/2019/07/sort-table-data-with-react/

var tableData = Object.values(run_data['sample_stats']);

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
  median_umis_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return a.Median_UMIs - b.Median_UMIs;
    }
  },
  median_umis_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return b.Median_UMIs - a.Median_UMIs;
    }
  },

  median_mito_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return a.Median_Mitochondrial_UMIs_Percent - b.Median_Mitochondrial_UMIs_Percents;
    }
  },
  median_mito_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return b.Median_Mitochondrial_UMIs_Percent - a.Median_Mitochondrial_UMIs_Percent;
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
  // doub_rate_up: {
  //   class: 'sort-up',
  //   fn: function fn(a, b) {
  //     return a.Doublet_Percent == "Fail" ? -1 : parseFloat(a.Doublet_Percent) - parseFloat(b.Doublet_Percent);
  //   }
  // },
  // doub_rate_down: {
  //   class: 'sort-down',
  //   fn: function fn(a, b) {
  //     return b.Doublet_Percent == "Fail" ? -1 : parseFloat(b.Doublet_Percent) - parseFloat(a.Doublet_Percent);
  //   }
  // },
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
  cfdr_p01_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return b.Cells_FDR_p01 - a.Cells_FDR_p01;
    }
  },
  cfdr_p01_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return a.Cells_FDR_p01 - b.Cells_FDR_p01;
    }
  },
  cfdr_p001_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return b.Cells_FDR_p001 - a.Cells_FDR_p001;
    }
  },
  cfdr_p001_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return a.Cells_FDR_p001 - b.Cells_FDR_p001;
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
    // }, _this.onSortDoub = function () {
    //   var currentSort = _this.state.currentSort;

    //   var nextSort = void 0;

    //   if (currentSort === 'doub_rate_down') nextSort = 'doub_rate_up';else if (currentSort === 'doub_rate_up') nextSort = 'doub_rate_down';else nextSort = 'doub_rate_up';

    //   _this.setState({
    //     currentSort: nextSort
    //   });
    }, _this.onSortUMIs = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'umis_down') nextSort = 'umis_up';else if (currentSort === 'umis_up') nextSort = 'umis_down';else nextSort = 'umis_up';

      _this.setState({
        currentSort: nextSort
      });

    }, _this.onSortMedianUMIs = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'median_umis_down') nextSort = 'median_umis_up';else if (currentSort === 'median_umis_up') nextSort = 'median_umis_down';else nextSort = 'median_umis_up';

      _this.setState({
        currentSort: nextSort
      });

    }, _this.onSortMedianMitoUMIs = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'median_mito_down') nextSort = 'median_mito_up';else if (currentSort === 'median_mito_up') nextSort = 'median_mito_down';else nextSort = 'median_umis_up';

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
    // }, _this.onSortc1000 = function () {
    //   var currentSort = _this.state.currentSort;

    //   var nextSort = void 0;

    //   if (currentSort === 'c1000_down') nextSort = 'c1000_up';else if (currentSort === 'c1000_up') nextSort = 'c1000_down';else nextSort = 'c1000_up';

    //   _this.setState({
    //     currentSort: nextSort
    //   });
    }, _this.onSortfdr_p01 = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'cfdr_p01_down') nextSort = 'cfdr_p01_up';else if (currentSort === 'cfdr_p01_up') nextSort = 'cfdr_p01_down';else nextSort = 'cfdr_p01_up';

      _this.setState({
        currentSort: nextSort
      });
    // }, _this.onSortfdr_p001 = function () {
    //   var currentSort = _this.state.currentSort;

    //   var nextSort = void 0;

    //   if (currentSort === 'cfdr_p001_down') nextSort = 'cfdr_p001_up';else if (currentSort === 'cfdr_p001_up') nextSort = 'cfdr_p001_down';else nextSort = 'cfdr_p001_up';

    //   _this.setState({
    //     currentSort: nextSort
    //   });
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
                "Median UMIs",
                React.createElement(
                  "button",
                  { onClick: this.onSortMedianUMIs, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Median Mitochondrial UMIs",
                React.createElement(
                  "button",
                  { onClick: this.onSortMedianMitoUMIs, className: "sort_button" },
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
              // React.createElement(
              //   "th",
              //   null,
              //   "Doublet rate",
              //   React.createElement(
              //     "button",
              //     { onClick: this.onSortDoub, className: "sort_button" },
              //     React.createElement("i", { className: "fas fa-sort" })
              //   )
              // ),
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
              // React.createElement(
              //   "th",
              //   null,
              //   "Cells with >1000 UMIs",
              //   React.createElement(
              //     "button",
              //     { onClick: this.onSortc1000, className: "sort_button" },
              //     React.createElement("i", { className: "fas fa-sort" })
              //   )
              // ),
              React.createElement(
                "th",
                null,
                "Cells with FDR<=.01",
                React.createElement(
                  "button",
                  { onClick: this.onSortfdr_p01, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              // React.createElement(
              //   "th",
              //   null,
              //   "Cells with FDR<=.001",
              //   React.createElement(
              //     "button",
              //     { onClick: this.onSortfdr_p001, className: "sort_button" },
              //     React.createElement("i", { className: "fas fa-sort" })
              //   )
              // )
            )
          ),
          React.createElement(
            "tbody",
            null,
            [].concat(_toConsumableArray(data)).sort(sortTypes[currentSort].fn).map(function (p, index) {
              return React.createElement(
                "tr",
                { key: index },
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
                  p.Median_UMIs
                ),
                React.createElement(
                  "td",
                  null,
                  p.Median_Mitochondrial_UMIs_Percent
                ),
                
                React.createElement(
                  "td",
                  null,
                  p.Duplication_rate
                ),
                // React.createElement(
                //   "td",
                //   null,
                //   p.Doublet_Percent
                // ),
                React.createElement(
                  "td",
                  null,
                  p.Cells_100_UMIs
                ),
                // React.createElement(
                //   "td",
                //   null,
                //   p.Cells_1000_UMIs
                // ),
                React.createElement(
                  "td",
                  null,
                  p.Cells_FDR_p01
                ),
                // React.createElement(
                //   "td",
                //   null,
                //   p.Cells_FDR_p001
                // )
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
  return React.createElement(Sample, { key: samp, sample_id: samp, garnett_model: run_data.sample_stats[samp].Garnett_model });
});

var sampPills = run_data.sample_list.map(function (samp) {
  return React.createElement(SamplePill, { key: samp, sample_id: samp });
});

ReactDOM.render(React.createElement(ExperimentPage, { samp_list: sampList, samp_pills: sampPills, run_name: run_data.run_name }), document.getElementById('exp_page'));
