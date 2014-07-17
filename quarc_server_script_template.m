%% Setup
% URI used to connect to the server model.
uri = 'shmem://foobar:1';
% Use blocking I/O. Do not change this value.
nonblocking = false;
% Connecting to the server using the specified URI
fprintf(1, 'Listening for clients at URI: %s\n', uri);
server = stream_listen(uri, nonblocking);
fprintf('Press Esc to exit this script. Do NOT press Ctrl+C.\n');
%% Establish connection
done = 0;
while ~done
    done = stream_poll(server, 0.1, 'accept');
    if qc_get_key_state(27) % Esc pressed stop
    end
end
%% Process control loop
try
    while ~qc_get_key_state(27) % while Esc key not pressed
        %% Receive data
        % Receive a double value from the client. This call will block
        % until data is received from the client or the client closes the
        % connection.
        %value = stream_receive_double(stream);
        if isempty(value) % then the client closed the connection gracefully
            fprintf(1, '\nClient has closed the connection.\n');
            break;
        end
        %% Send data
        % Store a double value in the stream send buffer
        %stream_send_double(stream, 2 * mod(t, 1));
        % Flush the send buffer to the underlying communications channel
        stream_flush(stream);
    end
    % Once the Esc key is pressed, close the stream handle used for
    % communications.
    fprintf(1, '\nShutting down the client connection...\n');
    stream_close(stream);
    fprintf(1, 'Client connection closed\n');
catch
    err = lasterror;
    fprintf(1, '\n%s.\nShutting down the client connection...\n', err.message);
    stream_close(stream);
    fprintf(1, 'Client connection closed\n');
    rethrow(err);
end
%% Close connection
% Once the Esc key is pressed, close the stream handle used for
% communications.
fprintf(1, '\nShutting down the server...\n');
stream_close(server);
fprintf(1, 'Server closed\n');