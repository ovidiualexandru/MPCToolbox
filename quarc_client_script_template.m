%% Setup
statedef;
% URI used to connect to the server model.
uri = 'shmem://foobar:1';
% Use blocking I/O. Do not change this value.
nonblocking = false;
% Connecting to the server using the specified URI
fprintf(1, 'Connecting to the server at URI: %s\n', uri);
stream = stream_connect(uri, nonblocking);
%Messages
fprintf(1, 'Connected to server.\n\n');
fprintf(1, 'Press Esc to exit this script. Do NOT press Ctrl+C.\n');
%% Process control loop
try
    while ~qc_get_key_state(27) % while Esc key not pressed
        %% Receive data
        % Receive a double value from the client. This call will block
        % until data is received from the client or the client closes the
        % connection.
        value = stream_receive_double_array(stream,8);
        %Do work
        fprintf('Received Data: Elevation: %f, Pitch %f\n', value(3), value(5));
        %% Send data
        % Store a double value in the stream send buffer
        u = struct;
        u.Vf = value(3);
        u.Vb = value(5);
        stream_send_array(stream, u);
        % Flush the send buffer to the underlying communications channel
        stream_flush(stream);
        if isempty(value) % then the server closed the connection gracefully
            fprintf(1, '\nServer has closed the connection.\n');
            break;
        end
    end
    % Once the Esc key is pressed, close the stream handle used for
    % communications.
    fprintf(1, '\nShutting down the client...\n');
    stream_close(stream);
    fprintf(1, 'Connection closed\n');
catch
    err = lasterror;
    fprintf(1, '\n%s.\nShutting down the client...\n', err.message);
    stream_close(stream);
    fprintf(1, 'Connection closed\n');
    rethrow(err);
end
%% Close connection
% Once the Esc key is pressed, close the stream handle used for
% communications.
fprintf(1, '\nShutting down the server...\n');
stream_close(server);
fprintf(1, 'Server closed\n');