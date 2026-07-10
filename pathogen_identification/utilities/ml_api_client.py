from typing import Any
from constants.constants import Televir_Metadata_Constants
import requests


class MLAPIClient:
    def __init__(self, timeout: int = 30):
        base_url = f"http://insaflu-ml-app:{Televir_Metadata_Constants.MODEL_PORT}"
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def _get(self, path: str, params: dict | None = None) -> dict[str, Any]:
        r = requests.get(f"{self.base_url}{path}", params=params, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def _post(self, path: str, body: dict | None = None) -> dict[str, Any]:
        r = requests.post(f"{self.base_url}{path}", json=body, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def health(self) -> dict[str, Any]:
        return self._get("/health")

    def models(self) -> dict[str, Any]:
        return self._get("/models")

    def reload(self, model_type: str | None = None) -> dict[str, Any]:
        if model_type:
            return self._post(f"/reload/{model_type}")
        return self._post("/reload")

    def predict_recall_cutoff(
        self,
        rows: list[dict[str, Any]],
        model: str = Televir_Metadata_Constants.RECALL_MODEL,
        tax_level: str = Televir_Metadata_Constants.RECALL_MODEL_TAX_LEVEL,
        target_recall: float | None = Televir_Metadata_Constants.TARGET_RECALL,
        confidence: float | None = None,
    ) -> dict[str, Any]:

        body: dict[str, Any] = {
            "model": model,
            "rows": rows,
            "tax_level": tax_level,
        }

        if target_recall is not None:
            body["target_recall"] = target_recall
        if confidence is not None:
            body["confidence"] = confidence

        return self._post("/predict_recall_cutoff_from_table", body)

    def predict_clustering_threshold(self, features: dict[str, Any]) -> dict[str, Any]:
        return self._post("/predict_televir_clustering_threshold", features)

    def predict_composition_stop_traversal(self, features: dict[str, float]) -> dict[str, Any]:
        return self._post("/predict_composition_stop_traversal", {"features": features})

