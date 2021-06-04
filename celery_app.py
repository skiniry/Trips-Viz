
from celery import Celery


celery_application = Celery("celery_test", backend='redis://127.0.0.1:6379', broker='redis://127.0.0.1:6379')
celery_application.autodiscover_tasks(['tripsviz.test_celery', 'tripsviz.find_orfs','tripsviz.generate_plot','tripsviz.generate_compare_plot'], force=True)

