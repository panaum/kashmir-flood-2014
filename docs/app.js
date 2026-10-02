/* Fills every number on the page from data/results.json, builds the map from
   data/layers.json, and runs the before/after slider. No figure is written in
   the HTML; no third-party request is made unless the visitor chooses the
   street map. */
(function () {
  "use strict";

  var MONTHS = ["January", "February", "March", "April", "May", "June", "July",
                "August", "September", "October", "November", "December"];

  function get(obj, path) {
    return path.split(".").reduce(function (o, k) { return o == null ? o : o[k]; }, obj);
  }

  function fmt(v, kind) {
    if (v == null) return "–";
    switch (kind) {
      case "int": return Math.round(v).toLocaleString("en-GB");
      case "dec1": return Number(v).toLocaleString("en-GB", { minimumFractionDigits: 1, maximumFractionDigits: 1 });
      case "dec2": return Number(v).toLocaleString("en-GB", { minimumFractionDigits: 2, maximumFractionDigits: 2 });
      case "dec3": return Number(v).toLocaleString("en-GB", { minimumFractionDigits: 3, maximumFractionDigits: 3 });
      case "date": {
        var m = String(v).match(/^(\d{4})-(\d{2})-(\d{2})/);
        return m ? (+m[3]) + " " + MONTHS[+m[2] - 1] + " " + m[1] : String(v);
      }
      default: return String(v);
    }
  }

  function fill(results) {
    document.querySelectorAll("[data-k]").forEach(function (el) {
      el.textContent = fmt(get(results, el.getAttribute("data-k")), el.getAttribute("data-fmt"));
    });
    document.querySelectorAll("[data-href]").forEach(function (el) {
      var href = get(results, el.getAttribute("data-href"));
      if (href) el.setAttribute("href", href);
    });
    var tbody = document.querySelector("#published tbody");
    results.published.figures.forEach(function (f) {
      var tr = document.createElement("tr");
      [f.source, f.figure.replace("km2", "km²"), f.date_range, f.sensor, f.difference].forEach(function (t, i) {
        var cell = document.createElement(i === 0 ? "th" : "td");
        if (i === 0) cell.setAttribute("scope", "row");
        cell.textContent = t;
        tr.appendChild(cell);
      });
      tbody.appendChild(tr);
    });
  }

  function setColours(layers) {
    var c = layers.layers.flood.colours, root = document.documentElement.style;
    root.setProperty("--c-flooded", c.flooded);
    root.setProperty("--c-persistent", c.persistent_water);
    root.setProperty("--c-lost", c.water_lost);
    root.setProperty("--c-na", c.not_assessed);
  }

  function buildMap(layers) {
    var L = window.L;
    var views = layers.views;
    var valley = L.latLngBounds(views.valley);
    var map = L.map("map", {
      zoomSnap: 0.25, keyboard: true, scrollWheelZoom: false,
      maxBounds: valley.pad(0.6), minZoom: 8, maxZoom: 16
    });
    map.attributionControl.setPrefix(false);

    var b = layers.layers.base, f = layers.layers.flood;
    var landsat = L.imageOverlay(b.file, b.bounds_latlon, {
      alt: "Landsat 8 false-colour image of the Kashmir valley, " + "25 August 2014",
      attribution: "Landsat 8 image courtesy of the U.S. Geological Survey"
    });
    var streets = null;
    var flood = L.imageOverlay(f.file, f.bounds_latlon, {
      opacity: 0.88, className: "flood-overlay",
      alt: "Flood classes: flooded, water on both dates, water lost, not assessed",
      attribution: "Flood map: this project; produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved"
    });
    landsat.addTo(map);
    flood.addTo(map);

    fetch(layers.layers.outline.file).then(function (r) { return r.json(); }).then(function (gj) {
      L.geoJSON(gj, { interactive: false, style: { color: "#F6F1E7", weight: 3.5, opacity: 0.8, fill: false } }).addTo(map);
      L.geoJSON(gj, { interactive: false, style: { color: "#2A1C13", weight: 1.5, dashArray: "5 4", fill: false } }).addTo(map);
    });

    map.fitBounds(views.srinagar);

    document.querySelectorAll(".views button").forEach(function (btn) {
      btn.addEventListener("click", function () {
        map.fitBounds(views[btn.getAttribute("data-view")]);
        document.querySelectorAll(".views button").forEach(function (o) {
          o.setAttribute("aria-pressed", String(o === btn));
        });
      });
    });

    document.querySelectorAll("input[name=base]").forEach(function (radio) {
      radio.addEventListener("change", function () {
        if (!radio.checked) return;
        if (radio.value === "osm") {
          if (!streets) {
            streets = L.tileLayer("https://tile.openstreetmap.org/{z}/{x}/{y}.png", {
              maxZoom: 19, referrerPolicy: "no-referrer",
              attribution: "Street map © <a href=\"https://www.openstreetmap.org/copyright\">OpenStreetMap</a> contributors"
            });
          }
          map.removeLayer(landsat);
          streets.addTo(map);
        } else {
          if (streets) map.removeLayer(streets);
          landsat.addTo(map);
          landsat.bringToBack();
        }
        flood.bringToFront();
      });
    });
  }

  function slider() {
    var range = document.getElementById("ba-range");
    var frame = document.querySelector(".ba-frame");
    if (!range || !frame) return;
    var set = function () {
      frame.style.setProperty("--pos", range.value + "%");
      range.setAttribute("aria-valuetext", range.value + " percent of the September image shown");
    };
    range.addEventListener("input", set);
    set();
  }

  slider();
  Promise.all([
    fetch("data/results.json").then(function (r) { return r.json(); }),
    fetch("data/layers.json").then(function (r) { return r.json(); })
  ]).then(function (d) {
    fill(d[0]);
    setColours(d[1]);
    buildMap(d[1]);
  }).catch(function () {
    var p = document.createElement("p");
    p.className = "note";
    p.textContent = "The data files for this page could not be loaded. If you opened the file directly, serve the docs/ folder over http instead.";
    document.getElementById("main").prepend(p);
  });
})();
