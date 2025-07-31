# app.py
from flask import Flask, render_template, request, jsonify
from core.predictor import KinaseMutationPredictor, load_kinase_list, load_gene_list, predict_once
import markdown

app = Flask(__name__)
predictor = KinaseMutationPredictor()
KINASES = load_kinase_list()
GENES   = load_gene_list()

@app.route("/")
def index():
    # Provide initial defaults similar to your Gradio examples
    defaults = {
        "kinase": KINASES[0] if KINASES else "",
        "gene": GENES[0] if GENES else "",
        "substrate": "RPQSPVGTGSY",  # 11-mer
        "position": 6,               # 1-based in UI
        "new_aa": "A"
    }
    return render_template("index.html", kinases=KINASES, genes=GENES, defaults=defaults)

@app.route("/api/predict", methods=["POST"])
def api_predict():
    data = request.get_json(force=True)
    kinase   = (data.get("kinase") or "").strip()
    gene     = (data.get("gene") or "").strip()
    substrate= (data.get("substrate") or "").strip().upper()
    position = int(data.get("position") or 0)  # 1-based
    new_aa   = (data.get("new_aa") or "").strip().upper()

    overall, details_md = predictor.predict_mutation_impact(
        kinase, gene, substrate, position-1, substrate[position-1] if substrate and 1 <= position <= len(substrate) else "", new_aa
    ) if substrate and 1 <= position <= len(substrate) else predict_once(kinase, gene, substrate, position, new_aa)

    details_html = markdown.markdown(details_md, extensions=["extra"])
    return jsonify({"overall": overall, "details_html": details_html})

if __name__ == "__main__":
    app.run(debug=True, host="127.0.0.1", port=5000)
