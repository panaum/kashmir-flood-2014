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

  /* ---- Atlas plate: cartouche, compass rose, scale bar, edge ticks ---- */

  var RESULTS = null;

  function atlasFurniture(map, L) {
    var Cartouche = L.Control.extend({
      options: { position: "topright" },
      onAdd: function () {
        var d = L.DomUtil.create("div", "cartouche");
        var date = RESULTS ? fmt(RESULTS.dates.flood, "date") : "";
        d.innerHTML = '<div class="ct-title">Kashmir valley</div>' +
          '<div class="ct-sub">Water on ' + date + '</div>' +
          '<div class="ct-sub">from Landsat 8</div>';
        return d;
      }
    });
    var Compass = L.Control.extend({
      options: { position: "bottomright" },
      onAdd: function () {
        var d = L.DomUtil.create("div", "compass");
        d.setAttribute("aria-hidden", "true");
        d.innerHTML = '<svg viewBox="0 0 64 64" fill="none" stroke="currentColor" stroke-width="1">' +
          '<circle cx="32" cy="34" r="20"/><circle cx="32" cy="34" r="16.5"/>' +
          '<path d="M32 12 L36 34 L32 56 L28 34 Z"/><path d="M10 34 L32 30 L54 34 L32 38 Z"/>' +
          '<path d="M32 12 L36 34 L28 34 Z" fill="currentColor"/>' +
          '<path d="M18 20 L33 33 M46 20 L31 33 M18 48 L31 35 M46 48 L33 35" stroke-width=".6"/>' +
          '<text x="32" y="8.5" text-anchor="middle" font-family="Bodoni Moda, serif" font-size="9" fill="currentColor" stroke="none">N</text>' +
          '</svg>';
        return d;
      }
    });
    new Cartouche().addTo(map);
    new Compass().addTo(map);
    var scaleMax = Math.max(70, Math.min(180, Math.round(map.getSize().x * 0.28)));
    L.control.scale({ position: "bottomleft", imperial: false, maxWidth: scaleMax }).addTo(map);

    var ticks = document.getElementById("ticks");
    var frame = ticks.parentNode;
    var STEPS = [60, 30, 20, 15, 10, 5, 2, 1, 0.5, 0.25];   // minutes of arc

    function label(v, isLat) {
      var hemi = isLat ? (v >= 0 ? "N" : "S") : (v >= 0 ? "E" : "W");
      var a = Math.abs(v), d = Math.floor(a + 1e-9), m = (a - d) * 60;
      var mm = Math.round(m * 100) / 100;
      if (mm >= 60) { d += 1; mm = 0; }
      var mins = Math.floor(mm + 1e-9), secs = Math.round((mm - mins) * 60);
      return d + "\u00b0" + (mins < 10 ? "0" : "") + mins + "\u2032" +
        (secs ? (secs < 10 ? "0" : "") + secs + "\u2033" : "") + "\u2009" + hemi;
    }

    function step(spanDeg, px) {
      var maxCount = Math.max(2, Math.floor(px / 110));
      for (var i = STEPS.length - 1; i >= 0; i--) {
        if (spanDeg * 60 / STEPS[i] <= maxCount) return STEPS[i] / 60;
      }
      return STEPS[0] / 60;
    }

    function el(cls, css, text) {
      var e = document.createElement("div");
      e.className = cls;
      Object.keys(css).forEach(function (k) { e.style[k] = css[k] + "px"; });
      if (text) e.textContent = text;
      ticks.appendChild(e);
      return e;
    }

    function draw() {
      ticks.innerHTML = "";
      var size = map.getSize(), b = map.getBounds(), c = b.getCenter();
      var mapEl = map.getContainer();
      var ox = mapEl.offsetLeft, oy = mapEl.offsetTop;     // map position inside the frame
      var len = 7;
      var dLon = step(b.getEast() - b.getWest(), size.x);
      for (var lon = Math.ceil(b.getWest() / dLon) * dLon; lon <= b.getEast(); lon += dLon) {
        var x = map.latLngToContainerPoint([c.lat, lon]).x;
        if (x < 24 || x > size.x - 24) continue;
        el("tick", { left: ox + x, top: oy - len, width: 1, height: len });
        el("tick", { left: ox + x, top: oy + size.y, width: 1, height: len });
        var t = el("tick-label", { left: ox + x, top: Math.max(2, oy - len - 14) }, label(lon, false));
        t.style.transform = "translateX(-50%)";
      }
      var dLat = step(b.getNorth() - b.getSouth(), size.y);
      for (var lat = Math.ceil(b.getSouth() / dLat) * dLat; lat <= b.getNorth(); lat += dLat) {
        var y = map.latLngToContainerPoint([lat, c.lng]).y;
        if (y < 24 || y > size.y - 24) continue;
        el("tick", { left: ox - len, top: oy + y, width: len, height: 1 });
        el("tick", { left: ox + size.x, top: oy + y, width: len, height: 1 });
        var s = el("tick-label", { left: ox - len - 2, top: oy + y }, label(lat, true));
        s.style.writingMode = "vertical-rl";
        s.style.transform = "translate(-100%, -50%) rotate(180deg)";
      }
    }
    map.on("move zoom resize viewreset load", draw);
    map.whenReady(draw);
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
      attribution: "Flood map: this project, using Copernicus WorldDEM-30 (<a href=\"#sources\">notice under Sources</a>)"
    });
    landsat.addTo(map);
    flood.addTo(map);

    fetch(layers.layers.outline.file).then(function (r) { return r.json(); }).then(function (gj) {
      L.geoJSON(gj, { interactive: false, style: { color: "#FFFFFF", weight: 3.5, opacity: 0.8, fill: false } }).addTo(map);
      L.geoJSON(gj, { interactive: false, style: { color: "#000000", weight: 1.5, dashArray: "5 4", fill: false } }).addTo(map);
    });

    atlasFurniture(map, L);
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
      range.setAttribute("aria-valuetext", range.value + " percent of the August image shown");
    };
    range.addEventListener("input", set);
    set();
  }

  slider();
  Promise.all([
    fetch("data/results.json").then(function (r) { return r.json(); }),
    fetch("data/layers.json").then(function (r) { return r.json(); })
  ]).then(function (d) {
    RESULTS = d[0];
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
