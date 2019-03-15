# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import os
import subprocess
import threading
import wx


class OcellarisClusterConnectorDialog(wx.Dialog):
    def __init__(self, parent, inspector_state, file_opener):
        super().__init__(
            parent,
            title='Cluster Connector',
            size=(550, 550),
            style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER,
        )
        self.istate = inspector_state
        self.file_opener = file_opener
        self.layout()

        self.avail_cluster_types = ['SLURM']
        self.handle_cluster = {'SLURM': (slurm_test, slurm_connect)}
        self.cluster_type.Set(self.avail_cluster_types)
        self.cluster_type.Select(0)

    def layout(self):
        p = wx.Panel(self)
        v = wx.BoxSizer(wx.VERTICAL)
        p.SetSizer(v)
        self.widgets1 = []
        self.widgets2 = []

        fg = wx.FlexGridSizer(cols=2, vgap=4, hgap=8)
        fg.AddGrowableCol(1, proportion=1)
        v.Add(fg, flag=wx.EXPAND | wx.ALL, border=5)

        fg.Add(wx.StaticText(p, label='Cluster type'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.cluster_type = wx.Choice(p)
        fg.Add(self.cluster_type, flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)

        fg.Add(wx.StaticText(p, label='Host name'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.host_name = wx.TextCtrl(p)
        fg.Add(self.host_name, flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)

        host_help = [
            'Either "user@host" or "host" depending on your setup.',
            'You MUST have password-less SSH login working!',
        ]
        fg.AddSpacer(5)
        fg.Add(wx.StaticText(p, label='\n'.join(host_help)))

        fg.Add(
            wx.StaticText(p, label='Verify connection'), flag=wx.ALIGN_CENTER_VERTICAL
        )
        b = wx.Button(p, label='Test')
        b.Bind(wx.EVT_BUTTON, self.test_connection)
        fg.Add(b, flag=wx.ALIGN_CENTER_VERTICAL)

        fg.AddSpacer(10)
        fg.Add(
            wx.StaticLine(p, style=wx.LI_HORIZONTAL),
            flag=wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL,
        )

        fg.Add(wx.StaticText(p, label='Cluster home'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.remote_dir = wx.TextCtrl(p)
        fg.Add(self.remote_dir, flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
        self.widgets1.append(self.remote_dir)

        fg.Add(wx.StaticText(p, label='Cluster user'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.remote_user = wx.TextCtrl(p)
        fg.Add(self.remote_user, flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
        self.widgets1.append(self.remote_user)

        fg.Add(wx.StaticText(p, label='Local mount'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.local_dir = wx.TextCtrl(p)
        fg.Add(self.local_dir, flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
        self.widgets1.append(self.local_dir)

        mount_help = [
            'The local directory where you have mounted the cluster',
            'home directory. E.g., "~/mymount".',
            'You can mount using "sshfs CLUSTER_HOST: ~/LOCAL_DIR"',
        ]
        fg.AddSpacer(5)
        fg.Add(wx.StaticText(p, label='\n'.join(mount_help)))

        fg.Add(
            wx.StaticText(p, label='Get running jobs'), flag=wx.ALIGN_CENTER_VERTICAL
        )
        b = wx.Button(p, label='Connect')
        b.Bind(wx.EVT_BUTTON, self.get_jobs)
        fg.Add(b, flag=wx.ALIGN_CENTER_VERTICAL)
        self.widgets1.append(b)

        fg.AddSpacer(10)
        fg.Add(
            wx.StaticLine(p, style=wx.LI_HORIZONTAL),
            flag=wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL,
        )

        # The file selector
        self.file_selector = wx.CheckListBox(p)
        v.Add(self.file_selector, flag=wx.EXPAND | wx.ALL, border=5, proportion=1)
        b = wx.Button(p, label='Load')
        b.Bind(wx.EVT_BUTTON, self.load_files)
        v.Add(b, flag=wx.ALIGN_RIGHT | wx.ALL, border=5)
        self.widgets2.append(self.file_selector)
        self.widgets2.append(b)

        for w in self.widgets1 + self.widgets2:
            w.Disable()

        p.Layout()

    def warning(self, title, message):
        with wx.MessageDialog(
            self, caption=title, message=message, style=wx.ICON_WARNING
        ) as dlg:
            dlg.ShowModal()

    def error(self, title, message):
        with wx.MessageDialog(
            self, caption=title, message=message, style=wx.ICON_EXCLAMATION
        ) as dlg:
            dlg.ShowModal()

    def test_connection(self, _evt=None):
        host = self.host_name.Value.strip()
        ic = self.cluster_type.GetSelection()
        cluster = self.avail_cluster_types[ic]
        tester = self.handle_cluster[cluster][0]

        homedir, username = tester(host)
        if homedir is None:
            self.error('Connection error', username)
            return

        for w in self.widgets1:
            w.Enable()

        self.remote_dir.Value = homedir
        self.remote_user.Value = username

    def get_jobs(self, _evt=None):
        host = self.host_name.Value.strip()
        user = self.remote_user.Value.strip()
        ic = self.cluster_type.GetSelection()
        cluster = self.avail_cluster_types[ic]
        connector = self.handle_cluster[cluster][1]

        directories = connector(host, user)
        if directories and directories[0] is None:
            self.error('Connection error', directories[1])
            return

        ok = self.process_directories(directories)
        if ok:
            self.file_selector.Set(self.files)

            for w in self.widgets2:
                w.Enable()

    def process_directories(self, directories):
        remote_dir = self.remote_dir.Value.strip()
        local_dir = self.local_dir.Value.strip()
        local_dir = os.path.expanduser(local_dir)
        N = len(remote_dir)

        self.files = []
        for d in directories:
            if not d.startswith(remote_dir):
                self.error(
                    'Path error',
                    'Directory %r\ndoes not start with\n%r' % (d, remote_dir),
                )
                return

            d2 = local_dir + d[N:]
            if not os.path.isdir(d2):
                self.error('Path error', 'Directory %r\ndoes not exist' % d2)
                return

            # Find the log file in this directory
            f = os.path.join(d2, 'ocellaris_out.log')
            if os.path.isfile(f):
                self.files.append(f)
            else:
                logs = []
                for f in os.listdir(d2):
                    if f.endswith('.log'):
                        logs.append(os.path.join(d2, f))
                if len(logs) > 1:
                    self.warning(
                        'Missing log file', 'Found multiple log files in\n' + d2
                    )
                elif len(logs) < 1:
                    self.warning('Missing log file', 'Found no log file in\n' + d2)
                else:
                    self.files.append(logs[0])

        return True

    def load_files(self, _evt=None):
        for i in self.file_selector.GetCheckedItems():
            f = self.files[i]
            self.file_opener(f)
        self.Close()


###############################################################################
# Cluster connectors


def slurm_test(host_name):
    cmd = ['ssh', host_name, 'pwd; whoami']
    retcode, stdout, _stderr, timed_out = run_command_with_timeout(cmd, 10)

    if timed_out:
        return None, 'Command timed out, did you enable passwordless login?'
    elif retcode != 0:
        return None, 'Command returned non-zero error code: %d' % retcode

    lines = stdout.strip().split('\n')
    if len(lines) < 2:
        return None, 'Did not get the expected output:\n%s' % stdout.strip()

    homedir = lines[-2]
    username = lines[-1]
    return homedir, username


def slurm_connect(host_name, user_name):
    cmd = ['ssh', host_name, 'squeue', '-o', '%Z', '-u', user_name]
    retcode, stdout, _stderr, timed_out = run_command_with_timeout(cmd, 10)

    if timed_out:
        return [None, 'Command timed out, did you enable passwordless login?']
    elif retcode != 0:
        return [None, 'Command returned non-zero error code: %d' % retcode]

    lines = stdout.strip().split('\n')
    return [line.strip() for line in lines[1:]]


def _kill_proc(proc, timeout):
    timeout["value"] = True
    proc.kill()


def run_command_with_timeout(cmd, timeout_sec):
    """
    Run a subprocess command with a timeout

    Returns (retcode, stdout, stderr, timed_out)

    https://stackoverflow.com/a/10768774/9293275
    """
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    timeout = {"value": False}
    timer = threading.Timer(timeout_sec, _kill_proc, [proc, timeout])
    timer.start()
    stdout, stderr = proc.communicate()
    timer.cancel()
    return (
        proc.returncode,
        stdout.decode("utf-8"),
        stderr.decode("utf-8"),
        timeout["value"],
    )
