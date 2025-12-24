"""
Logging configuration using loguru
"""
import sys
from loguru import logger
from app.config import settings

# Remove default handler
logger.remove()

# Disable all logging - only add file handler for critical errors
logger.add(
    "logs/app_{time:YYYY-MM-DD}.log",
    rotation="00:00",
    retention="30 days",
    level="CRITICAL",
    format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} - {message}",
)

__all__ = ["logger"]

