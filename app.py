import os

from flask import Flask
from flask_cors import CORS
from config import Config
from database import init_db

# Config configures MPLCONFIGDIR before Matplotlib is imported, avoiding a
# fallback to an unwritable user cache directory in worker processes.
import matplotlib
matplotlib.use('Agg')


_TRUE_VALUES = {'1', 'true', 'yes', 'on'}


def _env_enabled(name):
    """Return whether an opt-in runtime flag is enabled."""
    return os.environ.get(name, '').strip().lower() in _TRUE_VALUES


def _development_reload_files():
    """Return non-Python source files that should restart the dev server."""
    project_root = os.path.dirname(os.path.abspath(__file__))
    files = []
    for relative_root, extensions in (
        ('templates', {'.html'}),
        ('static', {'.css', '.js'}),
    ):
        source_root = os.path.join(project_root, relative_root)
        if not os.path.isdir(source_root):
            continue
        for directory, _, filenames in os.walk(source_root):
            files.extend(
                os.path.join(directory, filename)
                for filename in filenames
                if os.path.splitext(filename)[1].lower() in extensions
            )
    return files


def create_app():
    app = Flask(__name__)
    app.config.from_object(Config)

    from routes.auth import register_platform_access_gate
    register_platform_access_gate(app)

    cors_origins = os.environ.get('CORS_ORIGINS', '*')
    if cors_origins == '*':
        CORS(app)
    else:
        CORS(app, origins=cors_origins.split(','))

    init_db()

    from routes.main import main_bp
    from routes.analysis import analysis_bp
    from routes.projects import projects_bp
    from routes.upload import upload_bp
    from routes.results import results_bp
    from routes.figure_studio import figure_studio_bp
    from routes.api import api_bp
    from routes.chat import chat_bp
    from routes.branches import branches_bp
    from routes.workspace import workspace_bp
    from routes.ai_settings import ai_settings_bp
    from routes.wes import wes_bp

    app.register_blueprint(main_bp)
    app.register_blueprint(analysis_bp, url_prefix='/projects')
    app.register_blueprint(projects_bp, url_prefix='/projects')
    app.register_blueprint(upload_bp, url_prefix='/projects')
    app.register_blueprint(results_bp, url_prefix='/projects')
    app.register_blueprint(figure_studio_bp, url_prefix='/projects')
    app.register_blueprint(api_bp, url_prefix='/api')
    app.register_blueprint(chat_bp)
    app.register_blueprint(branches_bp)
    app.register_blueprint(workspace_bp, url_prefix='/projects')
    app.register_blueprint(ai_settings_bp)
    app.register_blueprint(wes_bp, url_prefix='/projects')

    return app

if __name__ == '__main__':
    app = create_app()
    port = int(os.environ.get('PORT', '5000'))
    auto_reload = _env_enabled('PLATFORM_AUTO_RELOAD')
    if auto_reload:
        # Keep the interactive debugger disabled on tunnel-exposed deployments.
        # The reloader is only a development convenience and restarts the whole
        # process, so it must not be enabled while long analyses are running.
        app.config['TEMPLATES_AUTO_RELOAD'] = True
        app.jinja_env.auto_reload = True
        print(
            'PLATFORM_AUTO_RELOAD 已启用：Python、模板和静态资源修改后将自动重载。',
            flush=True,
        )
    app.run(
        host='0.0.0.0',
        port=port,
        debug=False,
        threaded=True,
        use_reloader=auto_reload,
        extra_files=_development_reload_files() if auto_reload else None,
    )
