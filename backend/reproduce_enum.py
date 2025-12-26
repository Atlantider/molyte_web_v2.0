from enum import Enum
from typing import Dict, Any

class AccuracyLevel(str, Enum):
    STANDARD = "standard"

ACCURACY_CONFIGS = {
    AccuracyLevel.STANDARD: {"foo": "bar"}
}

def check(level):
    print(f"Checking {level} (type: {type(level)})")
    try:
        print(f"Result: {ACCURACY_CONFIGS[level]}")
    except KeyError as e:
        print(f"KeyError: {e}")

check(AccuracyLevel.STANDARD)
check("standard")
