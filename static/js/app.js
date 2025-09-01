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
      // Clear search inputs and restore all options
      el("kinase-search").value = '';
      el("gene-search").value = '';
      
      // Set values
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

// Search functionality for select dropdowns
function initSearchFilters() {
  // Store original options for both selects
  const kinaseSelect = el('kinase');
  const geneSelect = el('gene');
  const kinaseOptions = Array.from(kinaseSelect.options).map(opt => ({value: opt.value, text: opt.textContent}));
  const geneOptions = Array.from(geneSelect.options).map(opt => ({value: opt.value, text: opt.textContent}));

  function filterSelect(selectElement, options, searchTerm) {
    const currentValue = selectElement.value;
    selectElement.innerHTML = '';
    
    const filteredOptions = options.filter(option => 
      option.text.toLowerCase().includes(searchTerm.toLowerCase())
    );
    
    if (filteredOptions.length === 0) {
      const noResultOption = document.createElement('option');
      noResultOption.textContent = 'No matches found';
      noResultOption.disabled = true;
      selectElement.appendChild(noResultOption);
      return;
    }
    
    let hasCurrentValue = false;
    filteredOptions.forEach(option => {
      const optElement = document.createElement('option');
      optElement.value = option.value;
      optElement.textContent = option.text;
      if (option.value === currentValue) {
        optElement.selected = true;
        hasCurrentValue = true;
      }
      selectElement.appendChild(optElement);
    });
    
    // If current value is not in filtered results, select first option
    if (!hasCurrentValue && filteredOptions.length > 0) {
      selectElement.selectedIndex = 0;
      selectElement.dispatchEvent(new Event('change'));
    }
  }

  // Kinase search
  el('kinase-search').addEventListener('input', (e) => {
    filterSelect(kinaseSelect, kinaseOptions, e.target.value);
  });

  // Gene search  
  el('gene-search').addEventListener('input', (e) => {
    filterSelect(geneSelect, geneOptions, e.target.value);
  });

  // Add clear button functionality
  el('kinase-search').addEventListener('keydown', (e) => {
    if (e.key === 'Escape') {
      e.target.value = '';
      filterSelect(kinaseSelect, kinaseOptions, '');
    }
  });

  el('gene-search').addEventListener('keydown', (e) => {
    if (e.key === 'Escape') {
      e.target.value = '';
      filterSelect(geneSelect, geneOptions, '');
    }
  });

  // Clear search and restore all options when select changes
  kinaseSelect.addEventListener('change', () => {
    const searchInput = el('kinase-search');
    if (searchInput.value !== '') {
      searchInput.value = '';
      filterSelect(kinaseSelect, kinaseOptions, '');
    }
  });
  
  geneSelect.addEventListener('change', () => {
    const searchInput = el('gene-search');
    if (searchInput.value !== '') {
      searchInput.value = '';
      filterSelect(geneSelect, geneOptions, '');
    }
  });
}

["substrate","position","new_aa"].forEach(id=>{
  el(id).addEventListener("input", updatePreview);
  el(id).addEventListener("change", updatePreview);
});
el("predictBtn").addEventListener("click", predict);
bindExamples();
initSearchFilters();
updatePreview();
