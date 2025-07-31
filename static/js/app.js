const el = (id) => document.getElementById(id);

function updatePreview() {
  const substrate = el("substrate").value.toUpperCase();
  const pos = parseInt(el("position").value || "0", 10);
  const newAA = el("new_aa").value;
  const fromAA = (substrate && pos >= 1 && pos <= substrate.length) ? substrate[pos-1] : "?";
  const preview = el("preview");
  preview.querySelector(".from-aa").textContent = fromAA;
  preview.querySelector(".to-aa").textContent = newAA;
}

async function predict() {
  const btn = el("predictBtn");
  btn.classList.add("loading");
  try {
    const payload = {
      kinase: el("kinase").value,
      gene: el("gene").value,
      substrate: el("substrate").value.toUpperCase(),
      position: parseInt(el("position").value || "0", 10),
      new_aa: el("new_aa").value
    };
    const res = await fetch("/api/predict", {
      method: "POST",
      headers: {"Content-Type":"application/json"},
      body: JSON.stringify(payload)
    });
    const data = await res.json();
    el("overallImpact").textContent = data.overall || "—";
    el("detailsHtml").innerHTML = data.details_html || "";
    // expand details automatically on first result
    const acc = document.getElementById("detailsAcc");
    if (!acc.open) acc.open = true;
  } catch (e) {
    el("overallImpact").textContent = "Error";
    el("detailsHtml").innerHTML = "<p>Something went wrong. Check inputs and try again.</p>";
  } finally {
    btn.classList.remove("loading");
  }
}

function bindExamples() {
  document.querySelectorAll(".chip").forEach(chip => {
    chip.addEventListener("click", () => {
      el("kinase").value = chip.dataset.k;
      el("gene").value   = chip.dataset.g;
      el("substrate").value = chip.dataset.s;
      el("position").value  = chip.dataset.p;
      el("new_aa").value    = chip.dataset.aa;
      updatePreview();
      predict();
    });
  });
}

["substrate","position","new_aa"].forEach(id=>{
  el(id).addEventListener("input", updatePreview);
  el(id).addEventListener("change", updatePreview);
});
el("predictBtn").addEventListener("click", predict);
bindExamples();
updatePreview();
