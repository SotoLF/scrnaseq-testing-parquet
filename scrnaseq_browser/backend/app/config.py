import os, json
from pydantic import BaseModel
from typing import List, Dict, Any

class DatasetCfg(BaseModel):
    id: str
    name: str
    uri: str
    default_measurement: str = "RNA"
    default_obsm: str = "X_umap"

def get_datasets() -> List[DatasetCfg]:
    raw = os.getenv("DATASETS_JSON", "[]")
    data = json.loads(raw)
    return [DatasetCfg(**x) for x in data]

REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")
AWS_REGION = os.getenv("AWS_REGION", "us-east-2")
TILEDB_SOMA_THREADING = int(os.getenv("TILEDB_SOMA_THREADING", "1"))
