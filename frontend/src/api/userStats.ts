import axios from 'axios';

const API_BASE_URL = '/api/v1';

const getAuthHeaders = () => {
  const token = localStorage.getItem('access_token');
  return {
    Authorization: `Bearer ${token}`,
  };
};

export interface DailyStats {
  date: string;
  jobs_submitted: number;
  jobs_completed: number;
  jobs_failed: number;
  jobs_cancelled: number;
  cpu_hours_used: number;
  cluster_analysis_cpu_hours: number;
  cluster_analysis_task_count: number;
  max_concurrent_jobs: number;
  storage_used_gb: number;
}

/**
 * Get current user's daily usage statistics
 * @param days Number of days to retrieve (default: 7)
 */
export const getMyDailyStats = async (days: number = 7): Promise<DailyStats[]> => {
  const response = await axios.get(`${API_BASE_URL}/users/me/daily-stats`, {
    headers: getAuthHeaders(),
    params: { days },
  });
  return response.data;
};

/**
 * Get user's daily usage statistics (admin can view any user)
 * @param userId User ID
 * @param days Number of days to retrieve (default: 7)
 */
export const getUserDailyStats = async (userId: number, days: number = 7): Promise<DailyStats[]> => {
  const response = await axios.get(`${API_BASE_URL}/users/${userId}/daily-stats`, {
    headers: getAuthHeaders(),
    params: { days },
  });
  return response.data;
};

