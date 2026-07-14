from typing import Any
from constants.constants import Televir_Metadata_Constants
import requests


class MLAPIClient:
    def __init__(self, timeout: int = 30):
        base_url = f"http://insaflu-ml-app:{Televir_Metadata_Constants.MODEL_PORT}"
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def models_recall_cutoff(self) -> dict[str, Any]:
        """GET /models and filter"""
        models = self._get("/models").get('models', {})
        recall_cutoff_models = {}
        for model_info in models:
            if "recall" in model_info.get('model_type', ''):
                recall_cutoff_models[model_info.get('model_type', '')] = model_info
        return recall_cutoff_models
    
    def models_composition(self) -> dict[str, Any]:
        """GET /models and filter"""
        models = self._get("/models").get('models', {})
        composition_models = {}

        for model_info in models:
            if "composition" in model_info.get('model_type', ''):
                composition_models[model_info.get('model_type', '')] = model_info
        return composition_models

    def _get(self, path: str, params: dict | None = None) -> dict[str, Any]:
        r = requests.get(f"{self.base_url}{path}", params=params, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def _post(self, path: str, body: dict | None = None) -> dict[str, Any]:
        r = requests.post(f"{self.base_url}{path}", json=body, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def health(self) -> dict[str, Any]:
        """GET /health — liveness check."""
        return self._get("/health")

    def models(self) -> dict[str, Any]:
        """GET /models — list cached models with version and stage info."""
        return self._get("/models")

    def reload(self, model_type: str | None = None) -> dict[str, Any]:
        """POST /reload or /reload/{model_type} — reload model cache.

        Args:
            model_type: Composite key (e.g. ``"order_recall_gp_clf"``).
                        Omit to reload all models.
        """
        if model_type:
            return self._post(f"/reload/{model_type}")
        return self._post("/reload")

    def predict_recall_cutoff(
        self,
        rows: list[dict[str, Any]],
        model: str = "order_recall_gp_clf",
        target_recall: float | None = None,
        confidence: float | None = None,
    ) -> dict[str, Any]:
        """POST /predict_recall_cutoff_from_table — predict recall cutoff from raw table rows.

        Args:
            rows: List of dicts with keys ``taxid``, ``total_uniq_reads``, ``best_match_is_best``,
                  plus a column matching the model's ``tax_level`` (e.g. ``order``, ``family``).
            model: Full composite key from ``GET /models`` (e.g. ``"order_recall_gp_clf"``).
            target_recall: Optional target recall threshold.
            confidence: Optional confidence level for probability-guided cutoff.
        """
        body: dict[str, Any] = {
            "model": model,
            "rows": rows,
        }
        if target_recall is not None:
            body["target_recall"] = target_recall
        if confidence is not None:
            body["confidence"] = confidence
        return self._post("/predict_recall_cutoff_from_table", body)

    def predict_clustering_threshold(self, features: dict[str, Any]) -> dict[str, Any]:
        """POST /predict_televir_clustering_threshold — **deprecated**, returns 501."""
        return self._post("/predict_televir_clustering_threshold", features)

    def predict_composition_stop_traversal(self, features: dict[str, float], model: str = "xgb") -> dict[str, Any]:
        """POST /predict_composition_stop_traversal — predict stop_traversal from node features.

        Args:
            features: Dict of node feature values matching training column names.
            model: Composition variant (e.g. ``"xgb"``, ``"rf"``, ``"lr"``).
        """
        return self._post("/predict_composition_stop_traversal", {"model": model, "features": features})
