function Footer(props) {
  return React.createElement(
    "div",
    { className: "wrapper" },
    React.createElement(
      "div",
      { className: "footer-inner" },
      React.createElement(
        "div",
        { className: "footer-left" },
        React.createElement(
          "div",
          { className: "footer-logo" },
          React.createElement("img", { src: "https://brotmanbaty.org/wp-content/uploads/2018/06/BBI_Logo_Horizontal_Grey@2x.png", alt: "Footer Logo" })
        )
      ),
      React.createElement(
        "div",
        { className: "footer-right" },
        React.createElement(
          "div",
          { className: "fr-wrap" },
          React.createElement(
            "p",
            null,
            "357657 | Seattle, WA 98195-8047",
            React.createElement("br", null),
            "Health Sciences Building H-564",
            React.createElement("br", null),
            "info@brotmanbaty.org",
            React.createElement("br", null),
            "206-543-9660"
          )
        )
      )
    )
  );
}

function SubFooter(props) {
  return React.createElement(
    "div",
    { className: "wrapper" },
    React.createElement(
      "span",
      { className: "subfooter-logo sl-left" },
      React.createElement("img", { src: "https://brotmanbaty.org/wp-content/uploads/2018/06/UWMedLogoGrey@2x.png", alt: "Partner Logo" })
    ),
    React.createElement(
      "span",
      { className: "subfooter-logo sl-mid" },
      React.createElement("img", { src: "https://brotmanbaty.org/wp-content/uploads/2018/06/FredHutchLogoGrey@2x.png", alt: "Partner Logo" })
    ),
    React.createElement(
      "span",
      { className: "subfooter-logo sl-right" },
      React.createElement("img", { src: "https://brotmanbaty.org/wp-content/uploads/2018/06/SeattleChildrensLogoGrey@2x.png", alt: "Partner Logo" })
    )
  );
}

ReactDOM.render(React.createElement(Footer, null), document.getElementById('footer'));

ReactDOM.render(React.createElement(SubFooter, null), document.getElementById('subfooter'));