from __future__ import division
import os, json, re, threading, Queue, traceback
import cherrypy
from jinja2 import Environment, FileSystemLoader, select_autoescape

from plot_conc_space import parse_states, get_pfs, plot

FILE_DIR = os.path.abspath(os.path.dirname(__file__))
MEDIA_DIR = os.path.join(FILE_DIR, 'web')

jenv = Environment(
    loader=FileSystemLoader(os.path.join(MEDIA_DIR, 'templates')),
    autoescape=select_autoescape(['html', 'xml'])
)

templates = {
    'index': lambda: jenv.get_template('index.html'),
    'result': lambda: jenv.get_template('result.html'),
    'result_charts': lambda: jenv.get_template('result_charts.html'),
}

jobs_lock = threading.Lock()
job_file = lambda id: os.path.join(MEDIA_DIR, 'jobs', id + '.json')

try:
    os.mkdir(os.path.join(MEDIA_DIR, 'jobs'))
except:
    pass
    
class ChartPlugin(cherrypy.process.plugins.SimplePlugin):
    def plog(self, message):
        self.bus.log('[CHART PLUGIN]: {}'.format(message))
    
    def start(self):
        self.plog('Starting ChartPlugin')
        cherrypy.engine.subscribe('queue_job', self.queue_job)
        self.refresh_queue()
        
        self.plog('Starting plot generation')
        self.thread = threading.Thread(target=self.generate)
        self.thread.start()
    
    def refresh_queue(self):
        self.queue = Queue.Queue()
        self.plog('Refreshing queue')
        with jobs_lock:
            self.plog('Loading jobs file')
            try:
                with open(os.path.join(MEDIA_DIR, 'jobs.dat'), 'r') as f:
                    jobs = [{'id': job[0], 'email': job[1]} for job in [job.replace('\n', '').split(' ') for job in f.readlines()]]
            except IOError:
                jobs = []
        
        self.plog('Registering jobs')
        for job in jobs:
            self.queue.put(job)
    
    def stop(self):
        self.plog('Stopping ChartPlugin')
        # Remove all future tasks (will be restored from file on next boot)
        while not self.queue.empty():
            try:
                self.queue.get(False)
            except Empty:
                continue
            self.queue.task_done()
        # Wait until current task is done processing
        self.queue.put('STOP')
        self.queue.join()

    def queue_job(self, job_name, email):
        self.queue.put({'id': job_name, 'email': email })
            
    def generate(self):
        while True:
            job = self.queue.get(True)
            self.plog('job: {}'.format(job))
            if job == 'STOP':
                self.queue.task_done()
                return

            with open(job_file(job['id']), 'r') as f:
                job_def = json.load(f)
            self.plog('job def: {}'.format(job_def))

            try:
                for sequence in job_def['sequences']:
                    image_path = oss.path.join(MEDIA_DIR, 'static/plots/v1', job_def['round'], sequence + '.png')
                    if os.path.exists(image_path):
                        self.plog('{} exists'.format(image_path))
                        continue

                    try:
                        os.makedirs(os.path.dirname(image_path))
                    except OSError:
                        pass

                    self.plog('Parsing states')
                    inputs, reporter, complexes = parse_states(os.path.join(FILE_DIR, 'states', job_def['round'] + '.txt'))
                    self.plog('Running pfunc')
                    s = get_pfs(sequence, inputs, reporter, complexes)
                    self.plog('Plotting')
                    plot(s, image_path, nsteps=150, title=sequence)

                self.plog('Removing job')
                with jobs_lock:
                    with open(os.path.join(MEDIA_DIR, 'jobs.dat'), 'r+') as f:
                        lines = [line for line in f.readlines() if not line.startswith(job['id'])]
                        f.seek(0)
                        f.truncate()
                        f.write(''.join(lines))

                self.plog('Next')
                self.queue.task_done()
            except Exception, e:
                self.plog('Exception occurred, job will be skipped until next restart: {}'.format(e))
                with open(job_file(job['id']), 'r') as f:
                    job_def = json.load(f)
                job_def['error'] = traceback.format_exc()
                with open(job_file(job['id']), 'w') as f:
                    json.dump(job_def, f)
                self.queue.task_done()

ChartPlugin(cherrypy.engine).subscribe()

def get_result(job_id):
    with open(job_file(job_id)) as f:
        job = json.load(f)
    path = lambda sequence: os.path.join('v1', job['round'], sequence + '.png')
    charts = [{'src': path(sequence), 'sequence': sequence} for sequence in job['sequences'] if os.path.exists(os.path.join(MEDIA_DIR, 'static/plots', path(sequence)))]
    percent_complete = int(re.match('(\d+)', str((len(charts) / len(job['sequences']))*100)).group(0))
    
    return charts, percent_complete, job.get('error', None)

class App:
    def __init__(self):
        self.result = Result()
    
    @cherrypy.expose
    def index(self):
        return templates['index']().render(page='home')

    @cherrypy.expose
    def run_job(self, sequences, round, email=''):
        job_id = os.urandom(16).encode('hex')
        
        # It's a highly unlikely chance that this will run into a hash collision, but if we do, try again
        if os.path.exists(job_file(job_id)):
            return self.run_job(sequences, round, email)
        
        with open(job_file(job_id), 'w') as f:
            json.dump({'sequences': [seq.strip() for seq in sequences.split('\n') if seq.strip() != ''], 'round': round}, f)
        
        with jobs_lock:
            with open(os.path.join(MEDIA_DIR, 'jobs.dat'), 'a') as f:
                f.write('{} {}\n'.format(job_id, email))
                
        cherrypy.engine.publish('queue_job', job_id, email)
        
        raise cherrypy.HTTPRedirect('/result/' + job_id)
        
    @cherrypy.expose
    def jobstatus(self, job_id):
        charts, percent_complete, error = get_result(job_id)
        return json.dumps({'charts': templates['result_charts']().render(charts=charts), 'percent_complete': percent_complete, 'error': True if error is not None else False})

    @cherrypy.expose
    def dl_chart(self, path):
        return cherrypy.lib.static.serve_file(os.path.join(MEDIA_DIR, 'static/plots', path))
    
@cherrypy.popargs('job_id')
class Result:
    @cherrypy.expose
    def index(self, job_id):
        charts, percent_complete, error = get_result(job_id)
        return templates['result']().render(job=job_id, charts=charts, percent_complete=percent_complete, error=error)
    
if __name__ == '__main__':
    cherrypy.config.update({
        'server.socket_host': '0.0.0.0',
        'server.socket_port': 8081,
    })
    print('STATIC: {}'.format(os.path.join(MEDIA_DIR, 'static')))
    cherrypy.quickstart(App(), '/', config={
        '/static': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': os.path.join(MEDIA_DIR, 'static'),
        },
        '/favicon.ico':
        {
            'tools.staticfile.on': True,
            'tools.staticfile.filename': os.path.join(MEDIA_DIR, 'static/favicon.ico'),
        }
    })
