"""Local-only account administration. Run on the application server."""
import argparse
import getpass
import json
import os
import sqlite3
from pathlib import Path

from config import Config
from database import get_conn, init_db
from models import User
from werkzeug.security import generate_password_hash


def main():
    parser = argparse.ArgumentParser(description='平台用户管理')
    sub = parser.add_subparsers(dest='action', required=True)
    create = sub.add_parser('create')
    create.add_argument('username')
    create.add_argument('--admin', action='store_true')
    create.add_argument('--claim-existing', action='store_true', help='首次建管理员时认领所有历史项目')
    reset = sub.add_parser('reset-password')
    reset.add_argument('username')
    claim = sub.add_parser('claim-existing')
    claim.add_argument('username')
    for action in ('disable', 'enable'):
        item = sub.add_parser(action)
        item.add_argument('username')
    for action in ('grant-source-root', 'revoke-source-root'):
        item = sub.add_parser(action)
        item.add_argument('username')
        item.add_argument('assay', choices=('sc', 'wes'))
        item.add_argument('path')
    args = parser.parse_args()
    init_db()
    if args.action == 'create':
        if args.claim_existing and not args.admin:
            parser.error('--claim-existing 只能用于管理员账号')
        password = getpass.getpass('密码（至少 12 字符）：')
        if password != getpass.getpass('再次输入密码：'):
            parser.error('两次密码不一致')
        if len(password) < 12:
            parser.error('密码至少需要 12 个字符')
        if User.get_by_username(args.username):
            parser.error('用户名已存在')
        conn = get_conn()
        try:
            existing_admin = conn.execute('SELECT 1 FROM users WHERE is_admin=1 LIMIT 1').fetchone()
            unowned = conn.execute('SELECT COUNT(*) AS n FROM projects WHERE owner_user_id IS NULL').fetchone()['n']
        finally:
            conn.close()
        if args.claim_existing and existing_admin:
            parser.error('已有管理员；历史项目需要显式逐项核对归属')
        if unowned and not args.claim_existing:
            print(f'注意：{unowned} 个历史项目尚无归属，当前不可通过网页访问。')
        if args.claim_existing:
            # Keep a consistent copy of SQLite before changing ownership.
            db_file = Path(Config.DB_PATH)
            backup = db_file.with_suffix(db_file.suffix + '.before-user-ownership.bak')
            if backup.exists():
                parser.error(f'备份已存在：{backup}')
            source = sqlite3.connect(str(db_file))
            target = sqlite3.connect(str(backup))
            try:
                source.backup(target)
            finally:
                source.close()
                target.close()
            backup.chmod(0o600)
        user = User.create(args.username, password, is_admin=args.admin)
        if args.claim_existing:
            conn = get_conn()
            try:
                conn.execute('INSERT INTO user_ai_settings (user_id, key, value, updated_at) '
                             'SELECT ?, key, value, updated_at FROM platform_settings '
                             "WHERE key IN ('ai_api_url', 'ai_api_key', 'ai_model', 'ai_provider')", (user.id,))
                conn.execute("DELETE FROM platform_settings WHERE key IN ('ai_api_url', 'ai_api_key', 'ai_model', 'ai_provider')")
                conn.execute('UPDATE projects SET owner_user_id=? WHERE owner_user_id IS NULL', (user.id,))
                conn.commit()
            finally:
                conn.close()
            print(f'已认领 {unowned} 个历史项目；数据库备份：{backup}')
        print(f'已创建用户 {user.username}')
        return
    user = User.get_by_username(args.username)
    if not user:
        parser.error('用户不存在')
    if args.action == 'claim-existing':
        if not user.is_admin:
            parser.error('只能由管理员账号认领历史项目')
        db_file = Path(Config.DB_PATH)
        backup = db_file.with_suffix(db_file.suffix + '.before-user-ownership.bak')
        if backup.exists():
            parser.error(f'备份已存在：{backup}')
        source = sqlite3.connect(str(db_file))
        target = sqlite3.connect(str(backup))
        try:
            source.backup(target)
        finally:
            source.close()
            target.close()
        backup.chmod(0o600)
        conn = get_conn()
        try:
            count = conn.execute('SELECT COUNT(*) AS n FROM projects WHERE owner_user_id IS NULL').fetchone()['n']
            conn.execute('INSERT OR IGNORE INTO user_ai_settings (user_id, key, value, updated_at) '
                         'SELECT ?, key, value, updated_at FROM platform_settings '
                         "WHERE key IN ('ai_api_url', 'ai_api_key', 'ai_model', 'ai_provider')", (user.id,))
            conn.execute("DELETE FROM platform_settings WHERE key IN ('ai_api_url', 'ai_api_key', 'ai_model', 'ai_provider')")
            conn.execute('UPDATE projects SET owner_user_id=? WHERE owner_user_id IS NULL', (user.id,))
            conn.commit()
        finally:
            conn.close()
        print(f'已认领 {count} 个历史项目；数据库备份：{backup}')
        return
    if args.action in {'grant-source-root', 'revoke-source-root'}:
        path = os.path.realpath(args.path)
        if args.action == 'grant-source-root':
            if not os.path.isdir(path) or os.path.islink(args.path):
                parser.error('授权路径必须是现有非符号链接目录')
            roots = (Config.sc_batch_source_roots() if args.assay == 'sc' else
                     tuple(root for root in Config.wes_source_roots()
                           if root != os.path.realpath(Config.wes_upload_root())))
            if not any(path == root or path.startswith(root + os.sep) for root in roots):
                parser.error('授权目录必须位于管理员配置的源目录白名单内')
        column = 'sc_source_roots_json' if args.assay == 'sc' else 'wes_source_roots_json'
        grants = set(json.loads(getattr(user, column) or '[]'))
        if args.action == 'grant-source-root':
            grants.add(path)
        else:
            grants.discard(path)
        conn = get_conn()
        try:
            conn.execute(f'UPDATE users SET {column}=?, session_version=session_version+1 WHERE id=?',
                         (json.dumps(sorted(grants)), user.id))
            conn.commit()
        finally:
            conn.close()
        print(f'已更新 {user.username} 的 {args.assay} 源目录授权；该账号需重新登录')
        return
    conn = get_conn()
    try:
        if args.action == 'reset-password':
            password = getpass.getpass('新密码（至少 12 字符）：')
            if len(password) < 12:
                parser.error('密码至少需要 12 个字符')
            if password != getpass.getpass('再次输入新密码：'):
                parser.error('两次密码不一致')
            conn.execute('UPDATE users SET password_hash=?, session_version=session_version+1 WHERE id=?',
                         (generate_password_hash(password), user.id))
        else:
            conn.execute('UPDATE users SET is_active=?, session_version=session_version+1 WHERE id=?',
                         (int(args.action == 'enable'), user.id))
        conn.commit()
    finally:
        conn.close()
    print(f'已更新用户 {user.username}')


if __name__ == '__main__':
    main()
