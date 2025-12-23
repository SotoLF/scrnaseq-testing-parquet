import orjson
from typing import Optional, Any
import os
import zlib
from collections import OrderedDict
from threading import Lock

# Check if Redis should be used
USE_REDIS = os.getenv("USE_REDIS", "false").lower() == "true"

if USE_REDIS:
    import redis
    from .config import REDIS_URL
    r = redis.Redis.from_url(REDIS_URL, decode_responses=False)
else:
    # Local caching with guardrails:
    # - JSON cache for small objects
    # - Bytes cache (LRU) with strict size limits to avoid RAM spikes
    _json_cache = {}

    _bytes_cache_enabled = os.getenv("LOCAL_BYTES_CACHE", "true").lower() == "true"
    _bytes_cache_max_bytes = int(os.getenv("LOCAL_BYTES_CACHE_MAX_MB", "128")) * 1024 * 1024
    _bytes_cache_max_item_bytes = int(os.getenv("LOCAL_BYTES_CACHE_MAX_ITEM_MB", "16")) * 1024 * 1024

    _bytes_cache: "OrderedDict[str, bytes]" = OrderedDict()
    _bytes_cache_size = 0
    _lock = Lock()

def get_json(key: str) -> Optional[Any]:
    if USE_REDIS:
        val = r.get(key)
        if val is None:
            return None
        return orjson.loads(val)
    else:
        # Local cache
        if key in _json_cache:
            return _json_cache[key]
        return None

def set_json(key: str, obj: Any, ex: int = 3600) -> None:
    if USE_REDIS:
        r.set(key, orjson.dumps(obj), ex=ex)
    else:
        # Local cache (ignores expiration). Intentionally minimal.
        _json_cache[key] = obj

def get_bytes(key: str) -> Optional[bytes]:
    if USE_REDIS:
        return r.get(key)
    else:
        if not _bytes_cache_enabled:
            return None
        with _lock:
            val = _bytes_cache.get(key)
            if val is None:
                return None
            # LRU refresh
            _bytes_cache.move_to_end(key)
            return val

def set_bytes(key: str, b: bytes, ex: int = 3600) -> None:
    if USE_REDIS:
        r.set(key, b, ex=ex)
    else:
        if not _bytes_cache_enabled:
            return
        if b is None:
            return
        item_size = len(b)
        if item_size > _bytes_cache_max_item_bytes:
            return

        with _lock:
            global _bytes_cache_size
            # Replace existing
            if key in _bytes_cache:
                old = _bytes_cache.pop(key)
                _bytes_cache_size -= len(old)

            _bytes_cache[key] = b
            _bytes_cache.move_to_end(key)
            _bytes_cache_size += item_size

            # Evict LRU until within budget
            while _bytes_cache_size > _bytes_cache_max_bytes and len(_bytes_cache) > 0:
                _, old_val = _bytes_cache.popitem(last=False)
                _bytes_cache_size -= len(old_val)

def get_compressed(key: str) -> Optional[bytes]:
    """Get compressed bytes and decompress"""
    compressed = get_bytes(key)
    if compressed is None:
        return None
    try:
        return zlib.decompress(compressed)
    except Exception:
        return None

def set_compressed(key: str, data: bytes, ex: int = 3600) -> None:
    """Compress and store bytes. Use compression level 6 (good balance speed/size)"""
    compressed = zlib.compress(data, level=6)
    set_bytes(key, compressed, ex=ex)
