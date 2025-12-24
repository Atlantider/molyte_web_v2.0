"""
后端健康检查测试

验证后端 API 是否正常运行
"""

import pytest
from fastapi.testclient import TestClient
from app.main import app


@pytest.fixture
def client():
    """创建测试客户端"""
    return TestClient(app)


class TestBackendHealth:
    """后端健康检查测试"""

    def test_health_check(self, client):
        """测试健康检查端点"""
        response = client.get("/health")
        assert response.status_code == 200
        assert response.json() == {"status": "healthy"}

    def test_api_root(self, client):
        """测试 API 根路由"""
        response = client.get("/api/v1/")
        # 应该返回 404 或其他状态码，但不应该是 500
        assert response.status_code != 500

    def test_auth_me_without_token(self, client):
        """测试未授权的 /me 端点"""
        response = client.get("/api/v1/auth/me")
        # 应该返回 403 或 401
        assert response.status_code in [401, 403]

    def test_auth_me_with_invalid_token(self, client):
        """测试使用无效令牌的 /me 端点"""
        response = client.get(
            "/api/v1/auth/me",
            headers={"Authorization": "Bearer invalid_token"}
        )
        # 应该返回 401
        assert response.status_code == 401
        assert "detail" in response.json()


class TestOrganizationEndpoints:
    """组织 API 端点测试"""

    def test_organizations_endpoint_exists(self, client):
        """测试组织端点是否存在"""
        # 这个端点应该存在，但可能需要认证
        response = client.get("/api/v1/organizations")
        # 应该返回 401 或 200，但不应该是 404 或 500
        assert response.status_code != 404
        assert response.status_code != 500

    def test_organizations_admin_pending_endpoint_exists(self, client):
        """测试待审核组织端点是否存在"""
        response = client.get("/api/v1/organizations/admin/pending")
        # 应该返回 401 或 200，但不应该是 404 或 500
        assert response.status_code != 404
        assert response.status_code != 500


class TestCompensationEndpoints:
    """补偿 API 端点测试"""

    def test_compensation_endpoint_exists(self, client):
        """测试补偿端点是否存在"""
        response = client.get("/api/v1/compensation")
        # 应该返回 401 或 200，但不应该是 404 或 500
        assert response.status_code != 404
        assert response.status_code != 500


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

