import matplotlib
matplotlib.use('Agg')

from flask import Flask
from flask_cors import CORS
from config import Config
from database import init_db

def create_app():
    app = Flask(__name__)
    app.config.from_object(Config)
    CORS(app)

    init_db()

    from routes.main import main_bp
    from routes.projects import projects_bp
    from routes.upload import upload_bp
    from routes.analysis import analysis_bp
    from routes.results import results_bp
    from routes.api import api_bp

    app.register_blueprint(main_bp)
    app.register_blueprint(projects_bp, url_prefix='/projects')
    app.register_blueprint(upload_bp, url_prefix='/projects')
    app.register_blueprint(analysis_bp, url_prefix='/projects')
    app.register_blueprint(results_bp, url_prefix='/projects')
    app.register_blueprint(api_bp, url_prefix='/api')

    return app

if __name__ == '__main__':
    app = create_app()
    app.run(host='0.0.0.0', port=5000, debug=False, threaded=True)
