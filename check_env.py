import os
import json
from dotenv import load_dotenv

# Load .env manually to simulate what happens
load_dotenv()

raw = os.getenv("DATASETS_JSON", "[]")
print(f"Raw DATASETS_JSON: {raw}")

try:
    data = json.loads(raw)
    print(f"Parsed data: {json.dumps(data, indent=2)}")
except Exception as e:
    print(f"Error parsing JSON: {e}")
