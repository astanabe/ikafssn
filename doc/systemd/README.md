# systemd Unit Files for ikafssn

Sample systemd template unit files for running `ikafssnserver` and `ikafssnhttpd` as system services.

## Setup

### 1. Install unit files

```bash
sudo cp ikafssnserver@.service /etc/systemd/system/
sudo cp ikafssnhttpd@.service /etc/systemd/system/
sudo systemctl daemon-reload
```

### 2. Create configuration directory

```bash
sudo mkdir -p /etc/ikafssn
```

### 3. Create configuration files

For each database, create a server config and (optionally) an httpd config.

**Server config** (`/etc/ikafssn/nt.conf`):

```bash
IX_DIR=/data/ikafssn/nt_index
SOCKET_PATH=/run/ikafssn/nt.sock
THREADS=16
EXTRA_OPTS=
```

**HTTPD config** (`/etc/ikafssn/httpd-nt.conf`):

```bash
SERVER_SOCKET=/run/ikafssn/nt.sock
LISTEN=0.0.0.0:8080
THREADS=4
PATH_PREFIX=
EXTRA_OPTS=
```

### 4. Start services

```bash
# Start server
sudo systemctl enable --now ikafssnserver@nt

# Start HTTP frontend (requires server)
sudo systemctl enable --now ikafssnhttpd@nt
```

### 5. Check status

```bash
sudo systemctl status ikafssnserver@nt
sudo systemctl status ikafssnhttpd@nt
journalctl -u ikafssnserver@nt -f
```

## Multiple Databases

Use different instance names for each database:

```bash
# nt database
sudo cp /path/to/nt.conf /etc/ikafssn/nt.conf
sudo cp /path/to/httpd-nt.conf /etc/ikafssn/httpd-nt.conf
sudo systemctl enable --now ikafssnserver@nt ikafssnhttpd@nt

# refseq database
sudo cp /path/to/rs.conf /etc/ikafssn/rs.conf
sudo cp /path/to/httpd-rs.conf /etc/ikafssn/httpd-rs.conf
sudo systemctl enable --now ikafssnserver@rs ikafssnhttpd@rs
```

Use different listen ports or path prefixes for each httpd instance and configure nginx to route appropriately.

## TCP Mode

To use TCP instead of UNIX sockets, modify the server config:

```bash
# /etc/ikafssn/nt.conf
IX_DIR=/data/ikafssn/nt_index
SOCKET_PATH=
THREADS=16
EXTRA_OPTS=-tcp 0.0.0.0:9100
```

And the httpd config:

```bash
# /etc/ikafssn/httpd-nt.conf
SERVER_SOCKET=
LISTEN=0.0.0.0:8080
THREADS=4
PATH_PREFIX=
EXTRA_OPTS=-server_tcp 127.0.0.1:9100
```

Note: When using TCP mode, you need to adjust the `ExecStart` line in the unit file or override it with `systemctl edit`.
