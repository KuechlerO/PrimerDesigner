(() => {
  const KEY = "primerdesigner.uiScale";
  const MIN = 0.85;
  const MAX = 1.25;
  const STEP = 0.05;

  const clamp = (v) => Math.min(MAX, Math.max(MIN, v));
  const round2 = (v) => Math.round(v * 100) / 100;

  const getScale = () => {
    const raw = localStorage.getItem(KEY);
    const n = raw ? Number(raw) : NaN;
    return Number.isFinite(n) ? clamp(n) : 1;
  };

  const updateButtons = (scale) => {
    const down = document.getElementById("ui-scale-down");
    const up = document.getElementById("ui-scale-up");
    const reset = document.getElementById("ui-scale-reset");
    if (!down || !up || !reset) return;
    down.disabled = scale <= MIN + 1e-9;
    up.disabled = scale >= MAX - 1e-9;
    reset.disabled = Math.abs(scale - 1) < 1e-9;
  };

  const setScale = (v) => {
    const s = round2(clamp(v));
    document.documentElement.style.setProperty("--ui-scale", String(s));
    localStorage.setItem(KEY, String(s));
    updateButtons(s);
  };

  window.addEventListener("DOMContentLoaded", () => {
    // Apply early so layout uses the right em sizing.
    setScale(getScale());

    const down = document.getElementById("ui-scale-down");
    const up = document.getElementById("ui-scale-up");
    const reset = document.getElementById("ui-scale-reset");

    if (down) down.addEventListener("click", () => setScale(getScale() - STEP));
    if (up) up.addEventListener("click", () => setScale(getScale() + STEP));
    if (reset) reset.addEventListener("click", () => setScale(1));
  });
})();

