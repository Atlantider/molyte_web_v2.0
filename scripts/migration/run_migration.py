#!/usr/bin/env python3
import sys
import traceback
sys.path.insert(0, '/public/home/xiaoji/molyte_web/backend')

try:
    print("1. Importing database engine...")
    from app.database import engine
    from sqlalchemy import text
    print("✅ Database engine imported")

    # 读取迁移文件
    print("2. Reading migration file...")
    with open('/public/home/xiaoji/molyte_web/backend/migrations/add_resp_cpu_hours.sql', 'r') as f:
        sql = f.read()
    print(f"✅ Migration file read ({len(sql)} bytes)")

    # 执行迁移
    print("3. Executing migration...")
    with engine.connect() as conn:
        # 分割 SQL 语句并执行
        statements = [s.strip() for s in sql.split(';') if s.strip()]
        print(f"   Found {len(statements)} SQL statements")
        for i, stmt in enumerate(statements):
            try:
                conn.execute(text(stmt))
                print(f'   ✅ [{i+1}/{len(statements)}] {stmt[:60]}...')
            except Exception as e:
                print(f'   ⚠️  [{i+1}/{len(statements)}] {stmt[:60]}... - {str(e)[:100]}')
        conn.commit()
        print('✅ Migration completed!')
except Exception as e:
    print(f"❌ Error: {e}")
    traceback.print_exc()

