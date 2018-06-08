% Determines which computer the code is being run on - the iMac or the
% MacBook
function computer = comp_mac_address()

    localhost = java.net.InetAddress.getLocalHost;
    networkinterface = java.net.NetworkInterface.getByInetAddress(localhost);
    macaddress = typecast(networkinterface.getHardwareAddress, 'uint8');
    
    if min(macaddress == [56;201;134;2;215;75]) == 1
        computer = 'iMac';
    else, computer = 'MacBook';
    end

end