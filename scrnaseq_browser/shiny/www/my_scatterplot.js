function myScatterplotColors(payload, n) {
  const enc = payload && payload.colorEncoding;
  const c = payload && payload.color;

  if (!enc) {
    return new Array(n).fill('#4a90e2');
  }

  if (enc.type === 'categorical') {
    const levels = enc.levels || [];
    const numCategories = Math.max(1, levels.length);
    const palette = [];
    for (let i = 0; i < numCategories; i++) {
      const hue = (i * 360) / numCategories;
      palette.push(`hsl(${hue}, 70%, 60%)`);
    }
    const out = new Array(n);
    for (let i = 0; i < n; i++) {
      const code = c ? Number(c[i]) : 0;
      out[i] = palette[code] || '#999999';
    }
    return out;
  }

  if (enc.type === 'continuous') {
    const min = Number(enc.min);
    const max = Number(enc.max);
    const denom = (max - min) || 1;
    const out = new Array(n);
    for (let i = 0; i < n; i++) {
      const v = c ? Number(c[i]) : 0;
      const norm = Math.max(0, Math.min(1, (v - min) / denom));
      out[i] = viridisColor(norm);
    }
    return out;
  }

  return new Array(n).fill('#4a90e2');
}

function myScatterplotRender(el, payload, fallbackWidth, fallbackHeight) {
  if (!el || !payload) return;
  el.innerHTML = '';

  const xs = payload.x || [];
  const ys = payload.y || [];
  const n = Math.min(xs.length, ys.length);
  if (!n) return;

  const canvas = document.createElement('canvas');
  canvas.style.width = '100%';
  canvas.style.height = '100%';
  el.appendChild(canvas);

  const points = new Array(n);
  for (let i = 0; i < n; i++) {
    points[i] = [Number(xs[i]), Number(ys[i])];
  }

  const colors = myScatterplotColors(payload, n);
  const pointSize = Number(payload.pointSize || 3);
  renderScatterplot(canvas, el, points, colors, pointSize, fallbackWidth, fallbackHeight);
}

// Expose global renderer for plain Shiny rendering (no htmlwidgets)
window.my_scatterplot_render = myScatterplotRender;

// Optional: keep HTMLWidgets support if it's present
if (typeof HTMLWidgets !== 'undefined' && HTMLWidgets && typeof HTMLWidgets.widget === 'function') {
  HTMLWidgets.widget({
    name: 'my_scatterplot',
    type: 'output',
    factory: function(el, width, height) {
      return {
        renderValue: function(x) {
          // Accept both formats: old (x.points) and new (x.x/x.y)
          if (x && Array.isArray(x.points)) {
            const xs = x.points.map(p => Number(p.x));
            const ys = x.points.map(p => Number(p.y));
            const cc = x.points.map(p => p.color);
            myScatterplotRender(el, {
              x: xs,
              y: ys,
              color: cc,
              colorEncoding: x.colorEncoding,
              pointSize: x.pointSize
            }, width, height);
            return;
          }
          myScatterplotRender(el, x, width, height);
        },
        resize: function(width, height) {}
      };
    }
  });
}

// Simple canvas-based scatterplot renderer
function renderScatterplot(canvas, el, points, colors, pointSize, fallbackWidth, fallbackHeight) {
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  const { w, h, dpr } = (function() {
    const rect = el.getBoundingClientRect();
    const cw = rect.width || fallbackWidth || 800;
    const ch = rect.height || fallbackHeight || 400;
    const d = window.devicePixelRatio || 1;
    canvas.width = Math.max(1, Math.floor(cw * d));
    canvas.height = Math.max(1, Math.floor(ch * d));
    canvas.style.width = cw + 'px';
    canvas.style.height = ch + 'px';
    return { w: cw, h: ch, dpr: d };
  })();

  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  const bounds = (function() {
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (let i = 0; i < points.length; i++) {
      const x = points[i][0];
      const y = points[i][1];
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
    if (!isFinite(minX) || !isFinite(maxX) || !isFinite(minY) || !isFinite(maxY)) {
      minX = 0; maxX = 1; minY = 0; maxY = 1;
    }
    if (maxX === minX) maxX = minX + 1;
    if (maxY === minY) maxY = minY + 1;
    return { minX, maxX, minY, maxY };
  })();

  const padding = 20;
  const plotWidth = w - 2 * padding;
  const plotHeight = h - 2 * padding;

  // Clear canvas
  ctx.clearRect(0, 0, w, h);

  // Draw points
  points.forEach((point, i) => {
    const x = padding + ((point[0] - bounds.minX) / (bounds.maxX - bounds.minX)) * plotWidth;
    const y = padding + (1 - (point[1] - bounds.minY) / (bounds.maxY - bounds.minY)) * plotHeight;

    ctx.fillStyle = colors[i] || '#4a90e2';
    ctx.beginPath();
    ctx.arc(x, y, pointSize, 0, 2 * Math.PI);
    ctx.fill();
  });
}

// Viridis color scale approximation
function viridisColor(t) {
  // Simple viridis approximation
  const colors = [
    [68, 1, 84],
    [59, 82, 139],
    [33, 145, 140],
    [94, 201, 98],
    [253, 231, 37]
  ];

  const idx = Math.min(Math.floor(t * (colors.length - 1)), colors.length - 2);
  const localT = (t * (colors.length - 1)) - idx;

  const c1 = colors[idx];
  const c2 = colors[idx + 1];

  const r = Math.round(c1[0] + (c2[0] - c1[0]) * localT);
  const g = Math.round(c1[1] + (c2[1] - c1[1]) * localT);
  const b = Math.round(c1[2] + (c2[2] - c1[2]) * localT);

  return `rgb(${r}, ${g}, ${b})`;
}
