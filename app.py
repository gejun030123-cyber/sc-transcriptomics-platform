import matplotlib
matplotlib.use('Agg')

import os

from flask import Flask
from flask_cors import CORS
from config import Config
from database import init_db

def create_app():
    app = Flask(__name__)
    app.config.from_object(Config)

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
    from routes.api import api_bp
    from routes.chat import chat_bp

    app.register_blueprint(main_bp)
    app.register_blueprint(analysis_bp, url_prefix='/projects')
    app.register_blueprint(projects_bp, url_prefix='/projects')
    app.register_blueprint(upload_bp, url_prefix='/projects')
    app.register_blueprint(results_bp, url_prefix='/projects')
    app.register_blueprint(api_bp, url_prefix='/api')
    app.register_blueprint(chat_bp)

    return app

if __name__ == '__main__':
    app = create_app()
    app.run(host='0.0.0.0', port=5000, debug=False, threaded=True)
