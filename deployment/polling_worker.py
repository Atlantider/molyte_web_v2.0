#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Molyte Worker Entry Point (Refactored)
"""
import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

from worker.runner import CoreWorker

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Molyte Polling Worker")
    default_config = Path(__file__).parent / "polling_worker_config_tencent.yaml"
    parser.add_argument("--config", default=str(default_config), help="Path to config file")
    args = parser.parse_args()
    
    worker = CoreWorker(args.config)
    worker.run()
