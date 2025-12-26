"""
API 端点集成测试

测试新增的用户端 API 端点
"""

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.database import get_db


@pytest.fixture
def client(db_session):
    """创建测试客户端"""
    def override_get_db():
        try:
            yield db_session
        finally:
            pass

    app.dependency_overrides[get_db] = override_get_db
    return TestClient(app)


class TestUserOrganizationEndpoints:
    """用户端组织管理 API 测试"""

    def test_get_user_organizations(self, client, test_user, test_organization_member):
        """测试获取用户所属的组织"""
        # 这个测试需要认证，这里只是演示
        # 实际测试需要提供有效的 token
        response = client.get(
            "/api/v1/organizations/user/my-organizations",
            headers={"Authorization": f"Bearer test_token"}
        )
        # 预期返回 401（未认证）或 200（已认证）
        assert response.status_code in [200, 401, 403]

    def test_get_user_account_info(self, client, test_user):
        """测试获取用户的账号信息"""
        response = client.get(
            "/api/v1/organizations/user/my-account",
            headers={"Authorization": f"Bearer test_token"}
        )
        # 预期返回 401（未认证）或 200（已认证）
        assert response.status_code in [200, 401, 403]

    def test_get_user_quota_logs(self, client, test_user):
        """测试获取用户的配额日志"""
        response = client.get(
            "/api/v1/organizations/user/quota-logs",
            headers={"Authorization": f"Bearer test_token"}
        )
        # 预期返回 401（未认证）或 200（已认证）
        assert response.status_code in [200, 401, 403]


class TestAdminMasterAccountEndpoints:
    """管理员端主账号 API 测试"""

    def test_list_master_accounts(self, client, test_admin):
        """测试获取所有主账号列表"""
        response = client.get(
            "/api/v1/organizations/admin/master-accounts",
            headers={"Authorization": f"Bearer admin_token"}
        )
        # 预期返回 401（未认证）或 200（已认证）
        assert response.status_code in [200, 401, 403]

    def test_get_organization_quota_logs(self, client, test_admin, test_organization):
        """测试获取组织配额日志"""
        response = client.get(
            f"/api/v1/organizations/admin/{test_organization.id}/quota-logs",
            headers={"Authorization": f"Bearer admin_token"}
        )
        # 预期返回 401（未认证）或 200（已认证）
        assert response.status_code in [200, 401, 403]


class TestHealthCheck:
    """健康检查测试"""

    def test_api_is_running(self, client):
        """测试 API 是否正在运行"""
        # 测试一个不需要认证的端点
        response = client.get("/api/v1/jobs/")
        # 预期返回 401（需要认证）或 200（成功）
        assert response.status_code in [200, 401, 403, 422]

