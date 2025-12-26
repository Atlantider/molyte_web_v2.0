"""
API v1 router
"""
from fastapi import APIRouter
from app.api.v1 import auth, users, projects, electrolytes, jobs, slurm, admin, research, billing, visibility, qc, user_preferences, batch_import, worker, desolvation, binding, redox, cluster_analysis, forcefield, ai_discovery, similarity_search, compensation, accounts, recharge_packages, notifications, sms

api_router = APIRouter()

# Include all routers
api_router.include_router(auth.router, prefix="/auth", tags=["Authentication"])
api_router.include_router(users.router, prefix="/users", tags=["Users"])
api_router.include_router(accounts.router, tags=["Accounts"])
api_router.include_router(projects.router, prefix="/projects", tags=["Projects"])
api_router.include_router(electrolytes.router, prefix="/electrolytes", tags=["Electrolytes"])
api_router.include_router(jobs.router, prefix="/jobs", tags=["Jobs"])
api_router.include_router(qc.router, prefix="/qc", tags=["Quantum Chemistry"])
api_router.include_router(slurm.router, prefix="/slurm", tags=["Slurm"])
api_router.include_router(admin.router, tags=["Admin"])
api_router.include_router(research.router, prefix="/research", tags=["Research"])
api_router.include_router(billing.router, tags=["Billing"])
api_router.include_router(visibility.router, prefix="/visibility", tags=["Visibility"])
api_router.include_router(user_preferences.router, prefix="/user-preferences", tags=["User Preferences"])
api_router.include_router(batch_import.router, prefix="/batch-import", tags=["Batch Import"])
api_router.include_router(worker.router, tags=["Worker"])
api_router.include_router(desolvation.router, prefix="/desolvation", tags=["Desolvation Energy"])
api_router.include_router(binding.router, prefix="/binding", tags=["Binding Analysis"])
api_router.include_router(redox.router, tags=["Redox Potential"])
api_router.include_router(cluster_analysis.router, prefix="/cluster-analysis", tags=["Cluster Analysis"])
api_router.include_router(forcefield.router, prefix="/forcefield", tags=["Force Field"])
api_router.include_router(ai_discovery.router, prefix="/ai-discovery", tags=["AI Discovery"])
api_router.include_router(similarity_search.router, prefix="/similarity", tags=["Similarity Search"])
api_router.include_router(compensation.router, tags=["Compensation"])
api_router.include_router(recharge_packages.router, tags=["Recharge Packages"])
api_router.include_router(notifications.router, tags=["Notifications"])
api_router.include_router(sms.router, tags=["SMS"])
