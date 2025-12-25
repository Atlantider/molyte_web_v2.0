    def _update_job_cpu_hours_api(self, job_id: int, job_type: str, cpu_hours: float):
        """
        通过API更新任务的CPU hours
        
        Args:
            job_id: 任务ID
            job_type: 任务类型 ('qc', 'md', 'postprocess' 等)
            cpu_hours: CPU核时数
        """
        try:
            # 根据任务类型构建API路径
            type_mapping = {
                'qc': 'qc',
                'md': 'md',
                'postprocess': 'postprocess',
                'anion': 'anion-generation'
            }
            
            api_type = type_mapping.get(job_type, job_type)
            
            response = requests.patch(
                f"{self.api_base_url}/{api_type}/jobs/{job_id}",
                headers=self.api_headers,
                json={'actual_cpu_hours': cpu_hours},
                timeout=5
            )
            
            if response.status_code == 200:
                self.logger.debug(f"Updated {job_type} job {job_id} CPU hours: {cpu_hours:.2f}h")
            else:
                self.logger.debug(f"Failed to update CPU hours for {job_type} job {job_id}: {response.status_code}")
        except Exception as e:
            self.logger.debug(f"Error updating CPU hours for {job_type} job {job_id}: {e}")
