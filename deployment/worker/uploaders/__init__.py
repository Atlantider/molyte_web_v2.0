"""
Worker Uploaders 模块

包含COS上传、结果上传、部分结果上传
"""
from worker.uploaders.cos_client import COSClient
from worker.uploaders.result_uploader import ResultUploader
from worker.uploaders.partial_uploader import PartialUploader

__all__ = ['COSClient', 'ResultUploader', 'PartialUploader']
