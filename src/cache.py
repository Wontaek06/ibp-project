"""
Simple on-disk JSON cache for API responses.

Motivation: a full pipeline run makes ~100+ external API calls (GBIF,
Bio-ORACLE, UniProt, AlphaFold) and takes several minutes. During
development you re-run constantly while changing only the analysis code,
so caching the fetch layer makes iteration fast and avoids hammering
public APIs.

Usage:
    from src.cache import cached

    @cached("gbif")
    def gbif_records(species_name, limit=300):
        ...

Cache files land in .cache/<namespace>/<hash>.json. Delete the .cache
directory (or call clear_cache()) to force fresh fetches.

NOTE: the cache is keyed on the function arguments only. If you change
what a function *returns* (e.g. add a field), clear the cache or you'll
keep reading stale shapes.
"""
import functools
import hashlib
import json
import os
import shutil

CACHE_DIR = ".cache"


def _key(namespace, args, kwargs):
    payload = json.dumps({"a": args, "k": sorted(kwargs.items())}, default=str, sort_keys=True)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
    return os.path.join(CACHE_DIR, namespace, f"{digest}.json")


def cached(namespace, enabled=True):
    """Decorator: memoize a function's JSON-serializable return value to disk."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if not enabled or os.environ.get("IBP_NO_CACHE"):
                return func(*args, **kwargs)

            path = _key(namespace, args, kwargs)
            if os.path.exists(path):
                try:
                    with open(path, encoding="utf-8") as fh:
                        return json.load(fh)
                except (json.JSONDecodeError, OSError):
                    pass  # corrupt cache entry - just refetch

            result = func(*args, **kwargs)

            try:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w", encoding="utf-8") as fh:
                    json.dump(result, fh, ensure_ascii=False)
            except (TypeError, OSError):
                pass  # non-serializable or unwritable - skip caching silently

            return result
        return wrapper
    return decorator


def clear_cache(namespace=None):
    """Delete the whole cache, or one namespace."""
    target = os.path.join(CACHE_DIR, namespace) if namespace else CACHE_DIR
    if os.path.isdir(target):
        shutil.rmtree(target)
        print(f"Cleared cache: {target}")
    else:
        print(f"No cache at: {target}")
